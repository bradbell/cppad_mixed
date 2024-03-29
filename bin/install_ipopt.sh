#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-23 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ $0 != 'bin/install_ipopt.sh' ]
then
   echo 'bin/install_ipopt.sh: must be executed from its parent directory'
   exit 1
fi
#
# system_type
system_type="$1"
system_type_list=' debian red_hat mac_port mac_brew cygwin '
if ! echo "$system_type_list" | grep " $system_type " > /dev/null
then
   echo 'bin/install_ipopt.sh: system_type'
   echo "system_type='$system_type' is not one of the following:"
   echo "$system_type_list"
   exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# --------------------------------------------------------------------------
if [ "$system_type" == 'mac_brew' ]
then
   metis_libdir=$(brew --prefix)/lib
   metis_incdir=$(brew --prefix)/include
   if [ ! -e "$metis_libdir/libmetis.dylib" ]
   then
      echo 'mac_brew: Cannot find metis library directory'
   fi
   if [ ! -e "$metis_incdir/metis.h" ]
   then
      echo 'mac_brew: Cannot find metis include directory'
   fi
fi
# --------------------------------------------------------------------------
name_list=( \
   'ASL' \
   'Mumps' \
   'Ipopt' \
)
url_list=( \
   'https://github.com/coin-or-tools/ThirdParty-ASL' \
   'https://github.com/coin-or-tools/ThirdParty-Mumps' \
   'https://github.com/coin-or/Ipopt' \
)
version_list=( \
   'stable/2.0' \
   'stable/2.1' \
   'releases/3.13.4' \
)
java_flag_list=( \
   '' \
   '' \
   '--disable-java' \
)
# ---------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# specific_compiler
cmd=`grep '^specific_compiler=' bin/run_cmake.sh`
eval $cmd
#
# cmake_install_prefix
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
# --------------------------------------------------------------------------
# build_type.sh
bin/build_type.sh install_ipopt $build_type
# --------------------------------------------------------------------------
# cd external/build_type
if [ ! -e external/$build_type ]
then
   mkdir -p external/$build_type
fi
echo_eval cd external/$build_type
# -----------------------------------------------------------------------------
# n_job
if which nproc >& /dev/null
then
   n_job=$(nproc)
else
   n_job=$(sysctl -n hw.ncpu)
fi
# -----------------------------------------------------------------------------
# add_mumps_fcflags
#
# klugde necessary until coin or mumps fixes this problem
cat << EOF > junk.f
      program junk
      print*, "Hello World"
      end
EOF
if gfortran -c -fallow-argument-mismatch junk.f >& /dev/null
then
   echo 'Adding -fallow-argument-mismatch to Mumps fortran compiler flags'
   add_mumps_fcflags='ADD_FCFLAGS=-fallow-argument-mismatch'
else
   add_mumps_fcflags=''
fi
#-----------------------------------------------------------------------------
# debug_flags
if [ "$build_type" == 'debug' ]
then
   debug_flags='--enable-debug'
else
   debug_flags=''
fi
# -----------------------------------------------------------------------------
for i in {0..2}
do
   url=${url_list[$i]}
   version=${version_list[$i]}
   name=${name_list[$i]}
   java_flag=${java_flag_list[$i]}
   if [ ! -e $name.git ]
   then
      echo_eval git clone $url.git $name.git
   fi
   echo_eval cd $name.git
   echo_eval git reset --hard
   echo_eval git fetch origin
   echo_eval git checkout --quiet $version
   if [ -e "./get.$name" ]
   then
      if [ ! -e "./get.$name.done" ]
      then
         echo_eval ./get.$name
         touch ./get.$name.done
      fi
   fi
   if [ ! -e build ]
   then
      echo_eval mkdir build
   fi
   echo_eval cd build
   #
   add_fcflags=''
   with_metis_lflags=''
   with_metis_cflags=''
   metis_lflags=''
   metis_cflags=''
   if [ "$name" == 'Mumps' ]
   then
      add_fcflags="$add_mumps_fcflags"
      if [ "$system_type" == 'mac_brew' ]
      then
         with_metis_lflags='--with-metis-lflags='
         with_metis_cflags='--with-metis-cflags='
         metis_lflags="-L$metis_libdir -lmetis"
         metis_cflags="-I$metis_incdir"
      fi
   fi
cat << EOF
   ../configure \\
      --disable-dependency-tracking \\
      --prefix=$cmake_install_prefix \\
      --libdir=$cmake_install_prefix/$cmake_libdir \\
      --enable-shared \\
      $specific_compiler \\
      $java_flag $debug_flags $add_fcflags \\
      $with_metis_lflags"$metis_lflags" \\
      $with_metis_cflags"$metis_cflags"
EOF
   ../configure \
      --disable-dependency-tracking \
      --prefix=$cmake_install_prefix \
      --libdir=$cmake_install_prefix/$cmake_libdir \
      --enable-shared \
      $specific_compiler \
      $java_flag $debug_flags $add_fcflags \
      $with_metis_lflags"$metis_lflags" \
      $with_metis_cflags"$metis_cflags"
   echo_eval make -j $n_job install
   #
   # back to external
   echo_eval cd ../..
done
# ----------------------------------------------------------------------------
echo 'install_ipopt.sh: OK'
exit 0
