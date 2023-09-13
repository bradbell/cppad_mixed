#! /usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-23 Bradley M. Bell
# ----------------------------------------------------------------------------
set -e -u
#
if [ "$0" != "bin/build_type.sh" ]
then
   echo 'bin/build_type.sh: must be executed from its parent directory.'
   exit 1
fi
if [ "$#" != 2 ]
then
   cat << EOF
bin/build_type.sh program_name build_type

program_name:         is the program name for error reporting.
build_type:           this must be either debug or release
cmake_install_prefix: is its value in bin/run_cmake.sh

1. If uname -s begins with MSYS.* just warns that symbolic links are not set.
2. Otherwise sets up the following soft links:
   1. build -> build.built_type
   2. If cmake_install_prefix ends in /cppad_mixed:
      cmake_install_prefix -> cmake_install_prefix.build_type
EOF
   exit 1
fi
program_name="$1"
build_type="$2"
if [ "$build_type" != 'debug' ] && [ "$build_type" != 'release' ]
then
   echo "$program_name: build_type is not debug or release"
   exit 1
fi
# -----------------------------------------------------------------------------
kernel=$(uname -s)
if [[ "$kernel" =~ MSYS.* ]]
then
   echo 'Warning: MSYS does not suppor symbolic links'
   exit 0
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# -----------------------------------------------------------------------------
# cmake_install_prefix
eval $(grep '^cmake_install_prefix=' bin/run_cmake.sh)
#
# total_prefix
dollar='$'
total_prefix=`echo $cmake_install_prefix | sed \
   -e "s|/cppad_mixed$dollar|/cppad_mixed.$build_type|"`
echo "$total_prefix"
if [ "$total_prefix" == "$cmake_install_prefix" ]
then
   echo "$program_name: install prefix does not end in /cppad_mixed"
   echo 'Skipping soft link for install prefix'
else
# -----------------------------------------------------------------------------
   if [ ! -e "$total_prefix" ]
   then
      echo_eval mkdir "$total_prefix"
   fi
   if [ -e "$cmake_install_prefix" ]
   then
      if [ ! -L "$cmake_install_prefix" ]
      then
         echo "$program_name: $cmake_install_prefix is not a symbolic link"
         exit 1
      fi
      echo_eval rm "$cmake_install_prefix"
   fi
   echo_eval ln -s $total_prefix $cmake_install_prefix
# ---------------------------------------------------------------------------
fi
if [ -e build ]
then
   if [ ! -L build ]
   then
      echo "$program_name: build is not a symbolic link"
      exit 1
   fi
   echo_eval rm build
fi
if [ ! -e build.$build_type ]
then
   echo_eval mkdir build.$build_type
fi
echo_eval ln -s build.$build_type build
# ---------------------------------------------------------------------------
echo 'build_type.sh: OK'
exit 0
