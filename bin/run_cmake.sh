#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# =============================================================================
# $OMhelpKeyCharacter=&
# &begin run_cmake.sh&& &newlinech #&&
# &spell
#	cmake makefile cxx std dismod libdir lcppad cholmod xam cpp
#	usr eigen ipopt cppad bools suitesparse devel hpp Cholesky bool
# &&
#
# &section bin/run_cmake.sh: User Configuration Options&&
#
# &head cmake_verbose_makefile&&
# Use 'NO' for normal and 'YES' for verbose make output:
# &codep
cmake_verbose_makefile='NO'
# &&
#
# &head cmake_build_type&&
# Use either 'DEBUG' or 'RELEASE' for the type of this build:
# &codep
cmake_build_type='DEBUG'
# &&
#
# &head Prefixes&&
# Prefixes where the required packages are installed:
# &codep
cppad_prefix="$HOME/prefix/cppad_mixed"
eigen_prefix="$HOME/prefix/cppad_mixed/eigen"
ipopt_prefix="$HOME/prefix/cppad_mixed"
suitesparse_prefix="$HOME/prefix/cppad_mixed"
# &&
#
# &head extra_cxx_flags&&
# Extra C++ flags used during compilation:
# &codep
extra_cxx_flags='-std=c++11 -Wall'
# &&
#
# &head cmake_libdir&&
# Sub-directory of each prefix where libraries are installed.
# &codep
cmake_libdir='lib64'
# &&
#
# &head bool_sparsity&&
# If YES, use arrays of bools for sparsity patterns,
# otherwise use sets of indices.
# &codep
bool_sparsity='YES'
# &&
#
# &head cholmod_cholesky&&
# If YES, use &cref cholmod&& Cholesky factorization where possible.
# Otherwise always use Eigen's Cholesky factorization.
# &codep
cholmod_cholesky='YES'
# &&
#
# &head IHME Cluster Settings&&
# Here are some example changes that are used for the IHME cluster
# &codep
# extra_cxx_flags='-Wall'
# &&
# &end
# ============================================================================
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
if [ "$1" == '--help' ]
then
	cat << EOF
usage: bin/run_cmake.sh \\
	[--help] \\
	[--verbose] \\
	[--set_sparsity]
	[--eigen_cholesky]
EOF
	exit 0
else
	user_option="$1"
fi
if [ "$user_option" == '--verbose' ]
then
	cmake_verbose_makefile='1'
elif [ "$user_option" == '--set_sparsity' ]
then
	bool_sparsity='NO'
elif [ "$user_option" == '--eigen_cholesky' ]
then
	cholmod_cholesky='NO'
elif [ "$user_option" != '' ]
then
	echo "'$1' is an invalid option"
	bin/run_cmake.sh --help
	exit 1
fi
# ---------------------------------------------------------------------------
if [ ! -e build ]
then
	echo_eval mkdir build
fi
echo_eval cd build
if [ -e CMakeCache.txt ]
then
	rm CMakeCache.txt
fi
cmake \
	-Wno-dev \
	-D CMAKE_VERBOSE_MAKEFILE=$cmake_verbose_makefile \
	-D CMAKE_BUILD_TYPE=$cmake_build_type \
	\
	-D cppad_prefix="$cppad_prefix" \
	-D ipopt_prefix="$ipopt_prefix" \
	-D eigen_prefix="$eigen_prefix" \
	-D suitesparse_prefix="$suitesparse_prefix" \
	\
	-D extra_cxx_flags="$extra_cxx_flags" \
	-D cmake_libdir="$cmake_libdir" \
	-D bool_sparsity="$bool_sparsity" \
	-D cholmod_cholesky="$cholmod_cholesky" \
	..
# ---------------------------------------------------------------------------
echo 'run_cmake.sh: OK'
