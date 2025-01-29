#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ ! -e 'bin/check_install.sh' ]
then
   echo 'bin/check_install.sh: must be executed from its parent directory'
   exit 1
fi
#
# BEGIN_SORT_THIS_LINE_PLUS_3
# {xrst_begin check_install.sh}
# {xrst_spell
#     bitwise
#     cflags
#     cmd
#     config
#     cout
#     cp
#     endl
#     fi
#     grep
#     homebrew
#     int
#     lcppad
#     libdir
#     libs
#     mkdir
#     pkg
#     pkgconfig
#     rpath
#     sed
#     tmp
#     wl
#     wno
# }
# {xrst_comment END_SORT_THIS_LINE_MINUS_2}
# {xrst_comment_ch #}
#
# Example and Test Using the Installed Version of cppad_mixed
# ###########################################################
#
# Syntax
# ******
# ``bin/check_install.sh``
#
# build_type
# **********
# Get the :ref:`run_cmake.sh@build_type` used during the install:
# {xrst_code sh}
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
# {xrst_code}
#
# Prefixes
# ********
# Get the :ref:`prefixes<run_cmake.sh-name>` used during the install:
# {xrst_code sh}
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
ipopt_prefix="$cmake_install_prefix"
# {xrst_code}
#
# cmake_libdir
# ************
# Get the :ref:`run_cmake.sh@cmake_libdir`
# used during the install:
# {xrst_code sh}
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
# {xrst_code}
#
# example_file
# ************
# This is the user example that we will compile using
# the installed version of ``cppad_mixed`` :
# {xrst_code sh}
example_file='example/user/no_random.cpp'
# {xrst_code}
#
# PKG_CONFIG_PATH
# ***************
# Set the path used by $code pkg-config$$:
# {xrst_code sh}
if [ "$PKG_CONFIG_PATH" == '' ]
then
export PKG_CONFIG_PATH="$ipopt_prefix/$cmake_libdir/pkgconfig"
else
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$ipopt_prefix/$cmake_libdir/pkgconfig"
fi
# {xrst_code}
#
# LD_LIBRARY_PATH
# ***************
# This setting is no longer required so test with it empty:
# {xrst_code sh}
export LD_LIBRARY_PATH=''
# {xrst_code}
#
# Create Temporary
# ****************
# The following commands create a temporary directory,
# copy the example file into it,
# and make it the working directory:
# {xrst_code sh}
mkdir -p build/tmp
cp $example_file build/tmp/example.cpp
cd build/tmp
# {xrst_code}
#
# example_name
# ************
# The following command gets the example name, which is the example file
# with the ``.cpp`` at the end replaced by ``_xam`` :
# {xrst_code sh}
example_name=`echo $example_file | sed -e 's|.*/||' -e 's|\.cpp|_xam|'`
# {xrst_code}
#
# main
# ****
# The following command creates a main program
# (in the ``example.cpp`` file) that runs the example
# and reports the result of its return ``bool`` value:
# {xrst_code sh}
cat << EOF >> example.cpp
int main(void)
{  if( ! $example_name() )
   {  std::cout << "$example_name: Error" << std::endl;
      std::exit(1);
   }
   std::cout << "$example_name: OK" << std::endl;
   exit(0);
}
EOF
# {xrst_code}
#
# eigen_cflags
# ************
# The following command determines the flags necessary
# to compile with ``eigen`` on this system:
# {xrst_code sh}
eigen_cflags=$(pkg-config --cflags eigen3)
# {xrst_code}
#
# CHOLMOD_cflags
# **************
# The following command determines the flags necessary
# to compile with ``eigen`` on this system:
# {xrst_code sh}
CHOLMOD_cflags=$(pkg-config --cflags CHOLMOD)
# {xrst_code}
#
# gsl_libs
# ********
# The following command determines the library link flags necessary
# to link ``ipopt`` on this system:
# {xrst_code sh}
gsl_libs=`pkg-config --libs gsl`
# {xrst_code}
#
# ipopt_libs
# **********
# The following command determines the library link flags necessary
# to link ``ipopt`` on this system:
# {xrst_code sh}
ipopt_libs=`pkg-config --libs ipopt`
# {xrst_code}
#
# CHOLMOD_libs
# ************
# The following command sets the library link flags necessary
# to link ``cholmod`` on this system:
# {xrst_code sh}
CHOLMOD_libs=`pkg-config --libs CHOLMOD`
# {xrst_code}
#
# Compile and Link
# ****************
# The command below compiles and links the example program.
# {xrst_code sh}
#
# optimize_flags
if [ "$build_type" == 'debug' ]
then
   optimize_flags='-g -O0 -std=c++11 -Wall'
else
   optimize_flags='-O3 -DNDEBUG -std=c++11 -Wall'
fi
#
# path2libdir
path2libdir="$cmake_install_prefix/$cmake_libdir"
#
# homebrew_flags
if [ -e /opt/homebrew ]
then
   homebrew_flags='-Wno-bitwise-instead-of-logical'
   homebrew_flags+=' -I /opt/homebrew/include'
   homebrew_flags+=" -L /opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib"
else
   homebrew_flags=""
fi
#
cat << EOF
g++ example.cpp \\
   $optimize_flags \\
   -I $cmake_install_prefix/include \\
   -Wl,-rpath,$path2libdir -lcppad_mixed \\
   $homebrew_flags \\
   $eigen_cflags \\
   $CHOLMOD_cflags \\
   $gsl_libs \\
   $CHOLMOD_libs \\
   $ipopt_libs \\
   -o example
EOF
g++ example.cpp \
   $optimize_flags \
   -I $cmake_install_prefix/include \
   -Wl,-rpath,$path2libdir -lcppad_mixed \
   $eigen_cflags \
   $CHOLMOD_cflags \
   $homebrew_flags \
   $gsl_libs \
   $CHOLMOD_libs \
   $ipopt_libs \
   -o example
# {xrst_code}
#
# Run Example
# ***********
# The following commands run the example:
# {xrst_code sh}
if ! ./example
then
   echo 'check_install.sh: Error'
   exit 1
fi
echo 'check_install.sh: OK'
exit 0
# {xrst_code}
#
# {xrst_end check_install.sh}
