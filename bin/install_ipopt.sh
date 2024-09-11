#! /usr/bin/env bash
set -e -u
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-24 Bradley M. Bell
# ----------------------------------------------------------------------------
ipopt_url='https://github.com/coin-or/Ipopt'
ipopt_version='releases/3.13.4'
# ----------------------------------------------------------------------------
#
# echo_eval
function echo_eval() {
   echo $*
   eval $*
}
#
# clone_url_name_version
# creates ./name.git from corresponding url and version
# $1 = url
# $2 = name
# $3 = version
function clone_url_name_version() {
   if [ ! -e $2.git ]
   then
      echo_eval git clone $1.git $2.git
   fi
   echo_eval cd $2.git
   echo_eval git reset --hard
   echo_eval git fetch origin
   echo_eval git checkout --quiet $3
   if [ ! -e build ]
   then
      echo_eval mkdir build
   fi
   cd ..
}
# ----------------------------------------------------------------------------
if [ $0 != 'bin/install_ipopt.sh' ]
then
   echo 'bin/install_ipopt.sh: must be executed from its parent directory'
   exit 1
fi
#
# system_type
system_type="unknown"
if [ "$#" == 1 ]
then
   system_type="$1"
fi
system_type_list=' debian red_hat mac_port mac_brew cygwin '
if ! echo "$system_type_list" | grep " $system_type " > /dev/null
then
   echo 'bin/install_ipopt.sh: system_type'
   echo "system_type='$system_type' is not one of the following:"
   echo "$system_type_list"
   exit 1
fi
# ----------------------------------------------------------------------------
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
#
# build_type.sh
bin/build_type.sh install_ipopt $build_type
#
# external/$build_type
# cd external/build_type
if [ ! -e external/$build_type ]
then
   mkdir -p external/$build_type
fi
echo_eval cd external/$build_type
#
# n_job
if which nproc >& /dev/null
then
   n_job=$(nproc)
else
   n_job=$(sysctl -n hw.ncpu)
fi
#
# debug_flags
if [ "$build_type" == 'debug' ]
then
   debug_flags='--enable-debug'
else
   debug_flags=''
fi
#
# configure_all
configure_all="--disable-dependency-tracking"
configure_all+=" --prefix=$cmake_install_prefix"
configure_all+=" --libdir=$cmake_install_prefix/$cmake_libdir"
configure_all+=" --enable-shared"
configure_all+=" $specific_compiler"
configureall+=" $debug_flags"
#
# external/build_type/Ipopt.git
name='Ipopt'
clone_url_name_version $ipopt_url $name $ipopt_version
#
# external/build_type/ASL.git
# external/build_type/mumps.git
for name in 'ASL' 'Mumps'
do
   # clone_url_name_version
   line=$(grep "ThirdParty/$name" 'Ipopt.git/.coin-or/Dependencies')
   url=$(echo $line | awk '{print $2}' )
   version=$(echo $line | awk '{print $3}' )
   clone_url_name_version $url $name $version
   #
   # get.$name
   cd $name.git
   if [ -e "./get.$name" ]
   then
      if [ ! -e "./get.$name.done" ]
      then
         echo_eval ./get.$name
         touch ./get.$name.done
      fi
   fi
   cd ..
done
#
# Install ASL
cd ASL.git/build
echo_eval ../configure $configure_all
echo_eval make -j $n_job install
cd ../..
#
# Install Mumps
configure_mumps="$configure_all"
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
   configure_mumps+=" --with-metis-lflags='-L$metis_libdir -lmetis'"
   configure_mumps+=" --with-metis-cflags='-I$metis_incdir'"
fi
cd Mumps.git/build
echo_eval ../configure $configure_mumps
echo_eval make -j $n_job install
cd ../..
#
# Install Ipopt
cd Ipopt.git/build
echo_eval ../configure $configure_all --disable-java
echo_eval make -j $n_job install
cd ../..
# ----------------------------------------------------------------------------
echo 'install_ipopt.sh: OK'
exit 0
