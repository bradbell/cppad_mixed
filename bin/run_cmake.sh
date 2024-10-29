#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-24 Bradley M. Bell
# ----------------------------------------------------------------------------
#
# {xrst_begin run_cmake.sh}
# {xrst_spell
#     homebrew
#     wno
#     bitwise
#     callgrind
#     cd
#     config
#     cxx
#     fortran
#     gcc
#     gfortran
#     libdir
#     makefile
#     pc
#     pkg
#     wconversion
#     wpedantic
#     wshadow
# }
# {xrst_comment_ch #}
#
# bin/run_cmake.sh: User Configuration Options
# ############################################
#
# verbose_makefile
# ****************
# Use 'no' for normal and 'yes' for verbose make output:
# {xrst_code sh}
verbose_makefile='no'
# {xrst_code}
#
# build_type
# **********
# Use either 'debug' or 'release' for the type of this build:
# {xrst_code sh}
build_type='debug'
# {xrst_code}
#
# cmake_install_prefix
# ********************
# Prefix where cppad_mixed is installed:
# {xrst_code sh}
cmake_install_prefix="$HOME/prefix/cppad_mixed"
# {xrst_code}
# If *cmake_install_prefix* ends in ``/cppad_mixed`` ,
# ``run_cmake.sh`` will use a soft link from this prefix to
# *cmake_install_prefix* . ``debug`` or
# *cmake_install_prefix* . ``release``
# depending on the choice for *build_type* .
#
# Debug and Release
# *****************
# If a soft link is used for the install,
# the same technique will be used to map the ``build``
# directory to the debug or release version.
# If you are using both a debug and release versions of cppad_mixed,
# both versions of the
# :ref:`install_unix@Special Requirements`
# will need to be installed.
#
# Eigen Prefix
# ************
# It is a good idea to use a different prefix for installing eigen
# because the corresponding include files are treated like system files,
# otherwise the eigen include files would generate lots of warnings.
# The example install script ``bin/install_eigen.sh`` uses
# *cmake_install_prefix* / ``eigen`` as the prefix for installing eigen.
#
# specific_compiler
# *****************
# On some systems, e.g. the Mac using port, there are problems with mixing
# different compiler systems for fortran and C++; see
# `ipopt issue 471 <https://github.com/coin-or/Ipopt/discussions/471>`_.
# This variable allows you to set a specific compiler for
# C and or CXX and or FC.  For example
# ``specific_compiler='CC=gcc CXX=g++ FC=gfortran'``
# uses the gnu versions of these compilers.
# The configuration will automatically find compilers that are not specified;
# i.e., if
# {xrst_code sh}
specific_compiler=''
# {xrst_code}
#
# extra_cxx_flags
# ***************
# Extra C++ flags used to compile and test
# {xrst_code sh}
extra_cxx_flags='-Wpedantic -std=c++11 -Wall -Wshadow -Wconversion'
# for macOS using homebrew:
# 2DO: fix the warnings that are suppressed here.
# mac_brew: extra_cxx_flags+=' -I /opt/homebrew/include'
# mac_brew: extra_cxx_flags+=' -Wno-bitwise-instead-of-logical'
# mac_brew: extra_cxx_flags+=' -Wno-sign-conversion'
# {xrst_code}
#
# cmake_libdir
# ************
# Sub-directory of each prefix where libraries are installed.
# {xrst_code sh}
cmake_libdir='lib64'
# {xrst_code}
#
# cppad_mixed.pc
# ==============
# The file pkg-config file ``cppad_mixed.pc`` is installed in the
# *cmake_install_prefix/cmake_libdir* directory.
#
# ldlt_cholmod
# ************
# If yes, use ``ldlt_cholmod`` LDLT factorization where possible.
# Otherwise always use ``ldlt_eigen`` for LDLT factorization.
# {xrst_code sh}
ldlt_cholmod='yes'
# {xrst_code}
#
# optimize_cppad_function
# ***********************
# If yes, the operation sequence for certain CppAD functions
# will be optimized. This makes the code run faster but in some cases
# it can make debugging more complicated. It is suggested that you use
# ``no`` when *build_type* is ``debug`` and ``yes``
# when *build_type* is ``release`` .
# {xrst_code sh}
optimize_cppad_function='no'
# {xrst_code}
#
# for_hes_sparsity
# ****************
# If yes, user ``for_hes_sparsity`` to compute sparsity w.r.t. random
# effects (otherwise use ``rev_hes_sparsity`` ).
# {xrst_code sh}
for_hes_sparsity='yes'
# {xrst_code}
#
# Testing Speed and Memory
# ************************
# If you wish to test the speed or memory used by ``cppad_mixed`` ,
# set *build_type* to ``release`` , include the ``-g`` flag
# in *extra_cxx_flags* . Then execute the following commands:
#
# | |tab| ``bin/install_cppad.sh``
# | |tab| ``bin/run_cmake.sh``
# | |tab| ``cd build`` ; ``make`` *program* ; ``cd ..``
# | |tab| ``bin/`` *program* . ``sh`` *test2run*
#
# where *program* is
# :ref:`ar1_xam<ar1_xam.cpp-name>` or :ref:`capture_xam<ar1_xam.cpp-name>`
# and *test2run* is
# ``normal`` , ``callgrind`` , or ``massif`` .
#
# {xrst_end run_cmake.sh}
# ============================================================================
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# ----------------------------------------------------------------------------
if [ "$0" != "bin/run_cmake.sh" ]
then
   echo "bin/run_cmake.sh: must be executed from its parent directory"
   exit 1
fi
if [ "$build_type" != 'debug' ] && [ "$build_type" != 'release' ]
then
   echo "bin/run_cmake.sh: build_type is not 'debug' or 'release'"
   exit 1
fi
while [ "$1" != '' ]
do
   if [ "$1" == '--help' ]
   then
      cat << EOF
usage: bin/run_cmake.sh \\
   [--help] \\
   [--verbose] \\
   [--ldlt_eigen] \\
   [--rev_hes_sparsity] \\
   [--release] \\
   [--optimize_cppad_function] \\
   [--dismod_at_prefix]
EOF
      exit 0
   fi
   if [ "$1" == '--verbose' ]
   then
      verbose_makefile='1'
   elif [ "$1" == '--ldlt_eigen' ]
   then
      ldlt_cholmod='no'
   elif [ "$1" == '--rev_hes_sparsity' ]
   then
      for_hes_sparsity='no'
   elif [ "$1" == '--release' ]
   then
      build_type='release'
   elif [ "$1" == '--optimize_cppad_function' ]
   then
      optimize_cppad_function='yes'
   elif [ "$1" == '--dismod_at_prefix' ]
   then
      cmake_install_prefix="$HOME/prefix/dismod_at"
   else
      echo "'$1' is an invalid option"
      bin/run_cmake.sh --help
      exit 1
   fi
   shift
done
# Always set soft link for ./build -> ./build.buid_type
# If install prefix ends with cppad_mixed, also soft link install prefix
bin/build_type.sh run_cmake $build_type
# --------------------------------------------------------------------------
export PKG_CONFIG_PATH=':'
for pkg in ipopt cppad
do
   pkg_config_path=$(find -L "$cmake_install_prefix" -name "$pkg.pc" |\
      head -1 | sed -e "s|/$pkg.pc||" \
   )
   if [ "$pkg_config_path" == '' ]
   then
      echo "Can't find $pkg.pc: cmake_install_prefix=$cmake_install_prefix"
      exit 1
   fi
   if ! echo $PKG_CONFIG_PATH | grep ":$pkg_config_path:" > /dev/null
   then
      PKG_CONFIG_PATH=":$pkg_config_path$PKG_CONFIG_PATH"
   fi
done
PKG_CONFIG_PATH=$(echo $PKG_CONFIG_PATH | sed -e 's|^:||' -e 's|:$||')
echo "PKG_CONFIG_PATH=$PKG_CONFIG_PATH"
# --------------------------------------------------------------------------
# cmake_cxx_compiler
if echo $specific_compiler | grep 'CXX=' > /dev/null
then
   cxx=$(echo $specific_compiler | sed -e 's|.*CXX=\([^ ]*\).*|\1|')
   if ! which $cxx > /dev/null
   then
      echo "run_cmake.sh: specific_compiler: cannot execute $cxx compiler"
      exit 1
   fi
   cxx_path=$(which $cxx)
   cmake_cxx_compiler="-D CMAKE_CXX_COMPILER=$cxx_path"
else
   cmake_cxx_compiler=''
fi
# --------------------------------------------------------------------------
# cmake_c_compiler
if echo $specific_compiler | grep 'CC=' > /dev/null
then
   cc=$(echo $specific_compiler | sed -e 's|.*CC=\([^ ]*\).*|\1|')
   if ! which $cc > /dev/null
   then
      echo "run_cmake.sh: specific_compiler: cannot execute $cc compiler"
      exit 1
   fi
   c_path=$(which $cc)
   cmake_c_compiler="-D CMAKE_C_COMPILER=$c_path"
else
   cmake_c_compiler=''
fi
# --------------------------------------------------------------------------
if [ ! -e build ]
then
   echo_eval mkdir build
fi
echo_eval cd build
if [ -e CMakeCache.txt ]
then
   echo_eval rm CMakeCache.txt
fi
if [ -e CMakeFiles ]
then
   echo_eval rm -r CMakeFiles
fi
cat << EOF
cmake \\
   -Wno-dev \\
   -D CMAKE_VERBOSE_MAKEFILE=$verbose_makefile \\
   -D CMAKE_BUILD_TYPE=$build_type \\
   $cmake_cxx_compiler \\
   $cmake_c_compiler \\
   -D cmake_install_prefix="$cmake_install_prefix" \\
   \\
   -D extra_cxx_flags="$extra_cxx_flags" \\
   -D cmake_libdir="$cmake_libdir" \\
   -D ldlt_cholmod="$ldlt_cholmod" \\
   -D optimize_cppad_function="$optimize_cppad_function" \\
   -D for_hes_sparsity="$for_hes_sparsity" \\
   ..
EOF
cmake \
   -Wno-dev \
   -D CMAKE_VERBOSE_MAKEFILE=$verbose_makefile \
   -D CMAKE_BUILD_TYPE=$build_type \
   $cmake_cxx_compiler \
   \
   -D cmake_install_prefix="$cmake_install_prefix" \
   \
   -D extra_cxx_flags="$extra_cxx_flags" \
   -D cmake_libdir="$cmake_libdir" \
   -D ldlt_cholmod="$ldlt_cholmod" \
   -D optimize_cppad_function="$optimize_cppad_function" \
   -D for_hes_sparsity="$for_hes_sparsity" \
   ..
# ---------------------------------------------------------------------------
echo 'run_cmake.sh: OK'
exit 0
