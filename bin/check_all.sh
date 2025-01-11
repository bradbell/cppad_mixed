#! /usr/bin/env bash
set -e -u
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# -----------------------------------------------------------------------------
# echo_eval
echo_eval() {
   echo $*
   eval $*
}
#
# restore_and_exit
restore_and_exit() {
   $sed -i bin/run_cmake.sh \
      -e "s|^build_type=.*|build_type='debug'|" \
      -e 's|^\(.*\)\(# mac_brew:\) *|\2 \1|' \
      -e 's|  *$||'
   exit $1
}
# ----------------------------------------------------------------------------
if [ "$0" != "bin/check_all.sh" ]
then
   echo 'bin/check_all.sh must be executed from its parent directory'
   exit 1
fi
#
# grep, sed
source bin/grep_and_sed.sh
#
if ! $grep "build_type='debug'" bin/run_cmake.sh > /dev/null
then
   echo 'bin/check_all.sh: bin/run_cmake.sh build_type not debug'
   echo 'Attempted to automatically fix this.'
   restore_and_exit 1
fi
line_count=$($grep '^# mac_brew:' bin/run_cmake.sh | wc -l)
if [ $line_count != 3 ]
then
   echo 'check_all.sh: Cannot find three # mac_brew: lines in bin/run_cmake.sh'
   echo 'Attempted to automatically fix this.'
   restore_and_exit 1
fi
if ! $grep "optimize_cppad_function='no'" bin/run_cmake.sh > /dev/null
then
   echo 'bin/check_all.sh: bin/run_cmake.sh optimize_cppad_function not no'
   restore_and_exit 1
fi
#
if [ $# == 1 ]
then
   if [ "$1" == --help ]
   then
cat << EOF
bin/check_all.sh flags
possible flags
--ldlt_eigne               use eigen for LDLT factorization
--release                  only compile for release
--verbose_make             generate verbose makefiles
--skip_external_links      do not check documentation external links
--suppress_spell_warnings  do not check for documentaiton spelling errors
EOF
      exit 0
   fi
fi
#
# ldlt_eigen, release, verbose_make,
# skip_external_links, suppress_spell_warnings
ldlt_eigen='no'
release='no'
verbose_make='no'
skip_external_links='no'
suppress_spell_warnings='no'
while [ $# != 0 ]
do
   case "$1" in

      --ldlt_eigen)
      ldlt_eigen'yes'
      ;;

      --release)
      release='yes'
      ;;

      --verbose_make)
      verbose_make='yes'
      ;;

      --skip_external_links)
      skip_external_links='yes'
      ;;

      --suppress_spell_warnings)
      suppress_spell_warnings='yes'
      ;;

      *)
      echo "bin/check_all.sh: command line argument "$1" is not valid"
      exit 1
      ;;

   esac
   #
   shift
done
#
# bin/run_cmake.sh
if [ "$release" == 'yes' ]
then
   $sed -i bin/run_cmake.sh -e "s|^build_type=.*|build_type='release'|"
fi
if [ "$(uname)" == 'Darwin' ]
then
   if which brew
   then
      $sed -i -e 's|^\(# mac_brew:\) *\(.*\)|\2 \1|' bin/run_cmake.sh
   fi
fi
#
# bin/check_*.sh
list=`ls bin/check_*.sh`
for script in $list
do
   if [ "$script" != 'bin/check_all.sh' ] \
   && [ "$script" != 'bin/check_install.sh' ]
   then
      $script
   fi
done
#
# cmake_install_prefix
cmd=`$grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`$grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
if ls $cmake_install_prefix/$cmake_libdir/libcppad_mixed.* >& /dev/null
then
   echo_eval rm $cmake_install_prefix/$cmake_libdir/libcppad_mixed.*
fi
#
# installed_include_dir
installed_include_dir="$cmake_install_prefix/include/cppad/mixed"
if [ -e "$installed_include_dir" ]
then
   echo_eval rm -rf $installed_include_dir
fi
#
# run_cmake.sh
# The run_cmake.sh script sets up the build directory as a soft link.
# This should be done before xrst is run so it not create a normal directory
# called build.
flags=''
if [ "$release" == 'yes' ]
then
   flags="$flags --release"
fi
if [ "$ldlt_eigen" == 'yes' ]
then
   flags="$flags --ldlt_eigen"
fi
echo "bin/run_cmake.sh $flags >& cmake.log"
bin/run_cmake.sh $flags >& cmake.log
#
# run_xrst.sh
flags=''
if [ "$skip_external_links" == 'no' ]
then
   flags+=' --external_links'
fi
if [ "$suppress_spell_warnings" == 'yes' ]
then
   flags+=' --suppress_spell_warnings'
fi
bin/run_xrst.sh $flags
#
# build
cd build
#
# check
if which nproc >& /dev/null
then
   n_job=$(nproc)
else
   n_job=$(sysctl -n hw.ncpu)
fi
echo "make -j $n_job check >& check.log"
make -j $n_job check >& ../check.log
#
# speed
echo "make speed >& speed.log"
make speed >& ../speed.log
#
# install
echo "make -j $n_job install >& install.log"
echo_eval make -j $n_job install >& ../install.log
#
cd ..
#
if ! bin/check_install.sh
then
   echo 'bin/check_install.sh: Error'
   exit 1
fi
#
# warnings
# CppAD uses asserts to make sure this is not a problem
$sed -i check.log -e "/match_op.hpp:.*warning: ‘arg_match\[[01]\]’/d"
for target in cmake check speed install
do
   if $grep -i 'warning:' $target.log
   then
      echo "bin/run_check_all.sh: $target.log is has warnings."
      exit 1
   fi
done
# -----------------------------------------------------------------------------
echo 'check_all.sh: OK'
restore_and_exit 0
