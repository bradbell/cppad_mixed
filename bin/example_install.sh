#! /usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
set -e -u
#
#
# {xrst_begin example_install.sh}
# {xrst_comment_ch #}
#
# An Example Installation
# #######################
#
# Syntax
# ******
# ``bin/example_install.sh`` *run_test* *replace*
#
# run_test
# ********
# is either ``true`` or ``false`` .
# If it is true, this cppad_mixed tests will be run before installing.
# If there is an error in the tests, the install will abort.
#
# replace
# *******
# is either ``true`` or ``false`` .
# If an external is already installed and *replace* is true (false)
# the external will (will not) be replaced.
#
# Uninstall
# *********
# {xrst_toc_table
#     bin/example_remove.sh
# }
#
# Source
# ******
# {xrst_literal
#     BEGIN BASH
#     END BASH
# }
#
# {xrst_end example_install.sh}
#
# BEGIN BASH
if [ $0 != 'bin/example_install.sh' ]
then
   echo 'bin/example_install.sh: must be executed from its parent directory'
   exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# --------------------------------------------------------------------------
if [ "$#" != 2 ]
then
   echo 'bin/example_install.sh: run_test replace'
   echo 'run_test: is true or false'
   echo 'replace:  is true or false'
   exit 1
fi
if [ "$1" != 'true' ] && [ "$1" != 'false' ]
then
   echo 'bin/example_install.sh: run_test replace'
   echo 'run_test is not true or false'
   exit 1
fi
run_test="$1"
if [ "$2" != 'true' ] && [ "$2" != 'false' ]
then
   echo 'bin/example_install.sh: run_test replace'
   echo 'replace is not true or false'
   exit 1
fi
replace="$2"
# ---------------------------------------------------------------------------
# set build_type to value in run_cmake.sh
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
echo "build_type=$build_type"
#
# set cmake_install_prefix to value in run_cmake.sh
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
echo "cmake_install_prefix=$cmake_install_prefix"
#
# ipopt_prefix
ipopt_prefix="$cmake_install_prefix"
# --------------------------------------------------------------------------
# remove old version of example_install.log and example_install.err
for ext in log err
do
   if [ -e "example_install.$ext" ]
   then
      echo_eval rm "example_install.$ext"
   fi
done
# --------------------------------------------------------------------------
# set build link to build.debug or build.release depending on build_type
if echo "$cmake_install_prefix" | grep '/cppad_mixed$' > /dev/null
then
   bin/build_type.sh example_install.sh $build_type
fi
# --------------------------------------------------------------------------
user=$(whoami)
if [ "$user" == 'root' ]
then
   sudo=''
else
   sudo='sudo'
fi
# -----------------------------------------------------------------------------
# set system_type, system_install example_install.tmp
if which apt-get >& /dev/null
then
   system_type='debian'
   system_install="$sudo apt-get install -y"
   dpkg-query -l | sed -e 's|  *| |g' -e 's|^ii ||' > example_install.tmp
elif which dnf >& /dev/null
then
   system_type='red_hat'
   system_install="$sudo dnf install -y"
   dnf list --installed | sed -e 's|  *| |g' > example_install.tmp
elif which yum >& /dev/null
then
   system_type='red_hat'
   system_install="$sudo yum install -y"
   yum list installed | sed -e 's|  *| |g' > example_install.tmp
elif which port >& /dev/null
then
   system_type='mac_port'
   system_install="$sudo port install"
   port installed | sed -e 's|^ *||g' > example_install.tmp
elif which brew >& /dev/null
then
   system_type='mac_brew'
   system_install='brew install'
   brew list | sed -e 's|  *|\n|g' > example_install.tmp
elif which setup-x86_64 >& /dev/null
then
   system_tpye='cygwin'
   system_install='setup-x86_64.exe -q -P'
   cygcheck -c -d | sed -e 's|  *|-|' > example_install.tmp
else
   echo 'Cannot find the system pakcage manager'
   exit 1
fi
# --------------------------------------------------------------------------
# system external installs for normal system requirements
if [ "$system_type" == 'debian' ]
then
   list='
      cmake
      git
      wget
      libblas-dev
      libeigen3-dev
      liblapack-dev
      libsuitesparse-dev
      pkg-config
      g++
      gfortran
      libgsl-dev
   '
elif [ "$system_type" == 'red_hat' ]
then
   list='
      cmake
      eigen3-devel
      git
      wget
      blas-devel
      lapack-devel
      suitesparse-devel
      pkgconf
      gcc-c++
      gcc-gfortran
      gsl-devel
   '
elif [ "$system_type" == 'mac_port' ]
then
   list='
      cmake
      eigen3
      wget
      SuiteSparse
      pkgconfig
      gsl
   '
elif [ "$system_type" == 'mac_brew' ]
then
   # gnu-sed installs gsed
   # grep    installs ggrep
   list='
      gnu-sed
      grep
      cmake
      eigen
      wget
      suite-sparse
      pkg-config
      gsl
   '
elif [ "$system_tpye" == 'cygwin' ]
then
   list='
      gsl
      cmake
      eigen3
      git
      wget
      liblapack-devel
      pkgconf
      gcc-core
      gcc-g++
      gcc-fortran
      libgsl-devel
      patch
      libcholmod-devel
   '
else
   echo 'example_install.sh: script error'
   exit 1
fi
for package in $list
do
   if grep "^$package[^a-zA-Z_]" example_install.tmp > /dev/null
   then
      version=`grep "^$package[^a-zA-Z_]" example_install.tmp | head -1`
      echo "using installed $version"
   elif grep "^$package\$" example_install.tmp > /dev/null
   then
      # brew list case
      echo "using installed $package"
   else
      echo_eval $system_install $package
   fi
done
rm example_install.tmp
#
# ----------------------------------------------------------------------------
# local external installs for special requirements
for pkg in ipopt cppad
do
   # eval below converts $HOME in $prefix to its value for current user
   case $pkg in

      ipopt)
      eval file="$ipopt_prefix/include/coin-or/IpIpoptApplication.hpp"
      ;;

      cppad)
      eval file="$cmake_install_prefix/include/cppad/cppad.hpp"
      ;;

      *)
      echo 'bin/example_install.sh: program error'
      exit 1
      ;;
   esac
   #
   install='true'
   if [ -e "$file" ]
   then
      if [ "$replace" == 'false' ]
      then
         echo "Skipping bin/install_$pkg.sh"
         echo "Using previously installed version in $cmake_install_prefix"
         install='false'
      fi
   fi
   if [ "$install" == 'true' ]
   then
      name='example_install'
      echo "bin/install_$pkg.sh $system_type 1>> $name.log 2>> $name.err"
      if ! bin/install_$pkg.sh $system_type 1>> $name.log 2>> $name.err
      then
         echo "install_$pkg Error: Look at messages in"
         echo "tail $name.err"
         echo "tail $name.log"
         exit 1
      fi
   fi
done
# ----------------------------------------------------------------------------
# cppad_mixed
# ----------------------------------------------------------------------------
# bin/run_cmake.sh
echo "bin/run_cmake.sh 1>> example_install.log 2>> example_install.err"
if ! bin/run_cmake.sh 1>> example_install.log 2>> example_install.err
then
   tail example_install.err
   exit 1
fi
#
# change into build directory
echo_eval cd build
#
if which nproc >& /dev/null
then
   n_job=$(nproc)
else
   n_job=$(sysctl -n hw.ncpu)
fi
#
# make
if [ "$run_test" == 'true'  ]
then
   cmd_list='check speed install'
else
   cmd_list='install'
fi
for cmd in $cmd_list
do
   echo "make -j $n_job $cmd 1>> example_install.log 2>> example_install.err"
   if ! make -j $n_job $cmd \
      1>> ../example_install.log 2>> ../example_install.err
   then
      echo "Try following command in $(pwd) failed:"
      echo "    make -j $n_job $cmd"
      echo 'To see why look at example_install.err or run the comamnd'
      exit 1
   fi
done
cd ..
# ----------------------------------------------------------------------------
echo 'example_install.sh: OK'
exit 0
# END BASH
