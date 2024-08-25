#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-24 Bradley M. Bell
# -----------------------------------------------------------------------------
# echo_eval
echo_eval() {
   echo $*
   eval $*
}
#
# sed
if which gsed >& /dev/null
then
   sed=$(which gsed)
   echo "check_all.sh: replacing sed by $sed"
else
   sed='sed'
fi
# ----------------------------------------------------------------------------
if [ "$0" != "bin/check_all.sh" ]
then
   echo 'bin/check_all.sh must be executed from its parent directory'
   exit 1
fi
# ---------------------------------------------------------------------------
if ! grep "build_type='debug'" bin/run_cmake.sh > /dev/null
then
   echo 'bin/check_all.sh: bin/run_cmake.sh build_type not debug'
   exit 1
fi
if ! grep "optimize_cppad_function='no'" bin/run_cmake.sh > /dev/null
then
   echo 'bin/check_all.sh: bin/run_cmake.sh optimize_cppad_function not no'
   exit 1
fi
# ---------------------------------------------------------------------------
if [ "$1" == 'debug' ]
then
   release='no'
elif [ "$1" == 'release' ]
then
      release='yes'
      $sed -i bin/run_cmake.sh -e "s|^build_type=.*|build_type='release'|"
else
   echo 'usage: bin/check_all.sh (debug | release) [--ldlt_eigen]'
   exit 1
fi
if [ "$2" == '' ]
then
   ldlt_eigen='no'
elif [ "$2" == '--ldlt_eigen' ]
then
   ldlt_eigen='yes'
else
   echo 'usage: bin/check_all.sh (debug | release) [--ldlt_eigen]'
   exit 1
fi
# -----------------------------------------------------------------------------
#
# run_xrst.sh
bin/run_xrst.sh
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
# -----------------------------------------------------------------------------
# cmake_install_prefix
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
if ls $cmake_install_prefix/$cmake_libdir/libcppad_mixed.* >& /dev/null
then
   echo_eval rm $cmake_install_prefix/$cmake_libdir/libcppad_mixed.*
fi
#
installed_include_dir="$cmake_install_prefix/include/cppad/mixed"
if [ -e "$installed_include_dir" ]
then
   echo_eval rm -rf $installed_include_dir
fi
# -----------------------------------------------------------------------------
# run_cmake.sh
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
# ----------------------------------------------------------------------------
# build
# The run_cmake.sh script sets up the build directory as a soft link.
# This is necessary so xrst does not create a normal directory called build.
#
# build developer documentation
if [ -e build/html ]
then
   echo_eval rm -r build/html
fi
echo_eval xrst \
   --local_toc \
   --html_theme sphinx_rtd_theme \
   --group_list default dev
# ----------------------------------------------------------------------------
cd build
#
if which nproc >& /dev/null
then
   n_job=$(nproc)
else
   n_job=$(sysctl -n hw.ncpu)
fi
echo "make -j $n_job check >& check.log"
make -j $n_job check >& ../check.log
#
echo "make speed >& speed.log"
make speed >& ../speed.log
#
echo "make -j $n_job install >& install.log"
echo_eval make -j $n_job install >& ../install.log
#
cd ..
# -----------------------------------------------------------------------------
check_install_failed='no'
if [ "$release" == 'yes' ]
then
   if ! bin/check_install.sh
   then
      check_install_failed='yes'
   fi
   $sed -i bin/run_cmake.sh -e "s|^build_type=.*|build_type='debug'|"
else
   if ! bin/check_install.sh
   then
      check_install_failed='yes'
   fi
fi
if [ "$check_install_failed" == 'yes' ]
then
   echo 'bin/check_all.sh: Error'
   exit 1
fi
# -----------------------------------------------------------------------------
# CppAD uses asserts to make sure this is not a problem
$sed -i check.log -e "/match_op.hpp:.*warning: ‘arg_match\[[01]\]’/d"
for target in cmake check speed install
do
   if grep -i 'warning:' $target.log
   then
      echo "bin/run_check_all.sh: $target.log is has warnings."
      # 2DO: exit 1 after we fix some warning on MacOS
   fi
done
# -----------------------------------------------------------------------------
echo 'check_all.sh: OK'
exit 0
