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
#	usr eigen ipopt cppad bools suitesparse devel hpp ldlt bool
# &&
#
# &section bin/run_cmake.sh: User Configuration Options&&
#
# &head verbose_makefile&&
# Use 'NO' for normal and 'YES' for verbose make output:
# &codep
verbose_makefile='NO'
# &&
#
# &head build_type&&
# Use either 'debug' or 'release' for the type of this build:
# &codep
build_type='debug'
# &&
# Note that &code run_cmake.sh&& looks for a &icode%cppad_prefix%.debug%&&
# and &icode%cppad_prefix%.release%&& and uses them if they are present.
# Also note that it builds cppad_mixed in &code build.debug&& or
# &code build.release&& depending on &icode build_type&&.
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
# &head ldlt_cholmod&&
# If YES, use &cref ldlt_cholmod&& LDLT factorization where possible.
# Otherwise always use &cref ldlt_eigen&& for LDLT factorization.
# &codep
ldlt_cholmod='YES'
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
	[--release]
EOF
		exit 0
	fi
	if [ "$1" == '--verbose' ]
	then
		verbose_makefile='1'
	elif [ "$1" == '--ldlt_eigen' ]
	then
		ldlt_cholmod='NO'
	elif [ "$1" == '--release' ]
	then
		build_type='release'
	else
		echo "'$1' is an invalid option"
		bin/run_cmake.sh --help
		exit 1
	fi
	shift
done
# ---------------------------------------------------------------------------
if [ -e "$cppad_prefix.$build_type" ]
then
	if [ -e "$cppad_prefix" ]
	then
		echo_eval rm "$cppad_prefix"
	fi
	echo_eval ln -s $cppad_prefix.$build_type $cppad_prefix
fi
if [ ! -e "build.$build_type" ]
then
	echo_eval mkdir "build.$build_type"
fi
if [ -e build ]
then
	echo_eval rm build
fi
echo_eval ln -s build.$build_type build
echo_eval cd build
if [ -e CMakeCache.txt ]
then
	rm CMakeCache.txt
fi
cmake \
	-Wno-dev \
	-D CMAKE_VERBOSE_MAKEFILE=$verbose_makefile \
	-D CMAKE_BUILD_TYPE=$build_type \
	\
	-D cppad_prefix="$cppad_prefix" \
	-D ipopt_prefix="$ipopt_prefix" \
	-D eigen_prefix="$eigen_prefix" \
	-D suitesparse_prefix="$suitesparse_prefix" \
	\
	-D extra_cxx_flags="$extra_cxx_flags" \
	-D cmake_libdir="$cmake_libdir" \
	-D ldlt_cholmod="$ldlt_cholmod" \
	..
# ---------------------------------------------------------------------------
echo 'run_cmake.sh: OK'
