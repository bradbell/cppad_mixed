#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# =============================================================================
# $OMhelpKeyCharacter=&
# &begin run_cmake.sh&& &newlinech #&&
# &spell
#	cmake
#	makefile
#	cxx
#	std
#	dismod
#	libdir
#	lcppad
#	cholmod
#	cpp
#	usr
#	eigen
#	ipopt
#	cppad
#	bools
#	devel
#	hpp
#	ldlt
#	bool
#	libdir
#	cholesky
#	ar1_xam
#	callgrind
#	hes
#	Wshadow
#	Wconversion
#	Wpedantic
#	config
# &&
#
# &section bin/run_cmake.sh: User Configuration Options&&
#
# &head verbose_makefile&&
# Use 'no' for normal and 'yes' for verbose make output:
# &codep
verbose_makefile='no'
# &&
#
# &head build_type&&
# Use either 'debug' or 'release' for the type of this build:
# &codep
build_type='debug'
# &&
#
# &head cmake_install_prefix&&
# Prefix where cppad_mixed is installed:
# &codep
cmake_install_prefix="$HOME/prefix/cppad_mixed"
# &&
# If &icode cmake_install_prefix&& ends in &code /cppad_mixed&&,
# &code run_cmake.sh&& will use a soft link from this prefix to
# &icode%cmake_install_prefix%.debug%&& or
# &icode%cmake_install_prefix%.release%&&
# depending on the choice for &icode build_type&&.
#
# &head Debug and Release&&
# If a soft link is used for the install,
# the same technique will be used to map the &code build&&
# directory to the debug or release version.
# If you are using both a debug and release versions of cppad_mixed,
# both versions of the
# &cref/special requirements/install_unix/Special Requirements/&&
# will need to be installed.
#
# &head Eigen Prefix&&
# It is a good idea to use a different prefix for installing eigen
# because the corresponding include files are treated like system files,
# otherwise the eigen include files would generate lots of warnings.
# The example install script &code bin/install_eigen.sh&& uses
# &icode%cmake_install_prefix%/eigen%&& as the prefix for installing eigen.
#
# &head extra_cxx_flags&&
# Extra C++ flags used to compile and test
# &codep
extra_cxx_flags='-Wpedantic -std=c++11 -Wall -Wshadow -Wconversion'
# &&
#
# &head cmake_libdir&&
# Sub-directory of each prefix where libraries are installed.
# &codep
cmake_libdir='lib64'
# &&
#
# &subhead cppad_mixed.pc&&
# The file pkg-config file &code cppad_mixed.pc&& is installed in the
# &icode%cmake_install_prefix/cmake_libdir%&& directory.
#
# &head ldlt_cholmod&&
# If yes, use &code ldlt_cholmod&& LDLT factorization where possible.
# Otherwise always use &code ldlt_eigen&& for LDLT factorization.
# &codep
ldlt_cholmod='yes'
# &&
#
# &head optimize_cppad_function&&
# If yes, the operation sequence for certain CppAD functions
# will be optimized. This makes the code run faster but in some cases
# it can make debugging more complicated. It is suggested that you use
# &code no&& when &icode build_type&& is &code debug&& and &code yes&&
# when &icode build_type&& is &code release&&.
# &codep
optimize_cppad_function='no'
# &&
#
# &head for_hes_sparsity&&
# If yes, user &code for_hes_sparsity&& to compute sparsity w.r.t. random
# effects (otherwise use &code rev_hes_sparsity&&).
# &codep
for_hes_sparsity='yes'
# &&
#
# &head Testing Speed and Memory&&
# If you wish to test the speed or memory used by &code cppad_mixed&&,
# set &icode build_type&& to &code release&&, include the &code -g&& flag
# in &icode extra_cxx_flags&&. Then execute the following commands:
# &codei%
#	bin/install_cppad.sh
#	bin/run_cmake.sh
#	cd build; make %program%; cd ..
#	bin/%program%.sh %test2run%
# %&&
# where &icode program&& is
# &cref/ar1_xam/ar1_xam.cpp/&& or &cref/capture_xam/ar1_xam.cpp/&&
# and &icode test2run&& is
# &code normal&&, &code callgrind&&, or &code massif&&.
#
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
# If install prefix end is cppad_mixed, also set soft link for install dir
bin/build_type.sh run_cmake $cmake_install_prefix $build_type
# --------------------------------------------------------------------------
export PKG_CONFIG_PATH=':'
for pkg in eigen3 ipopt cppad
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
cmake \
	-Wno-dev \
	-D CMAKE_VERBOSE_MAKEFILE=$verbose_makefile \
	-D CMAKE_BUILD_TYPE=$build_type \
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
