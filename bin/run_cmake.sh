#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
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
#	suitesparse
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
# Note that if &icode cppad_prefix&& ends in &code cppad_mixed&&,
# &code run_cmake.sh&& will use a link from the prefix to
# &icode%cppad_prefix%.debug%&& or
# &icode%cppad_prefix%.release%&&
# depending on the choice &icode build_type&&.
#
# &head Prefixes&&
# Prefixes where the required packages are installed:
# &codep
cppad_prefix="$HOME/prefix/cppad_mixed"
eigen_prefix="$HOME/prefix/cppad_mixed/eigen"
ipopt_prefix="$HOME/prefix/cppad_mixed"
suitesparse_prefix="$HOME/prefix/cppad_mixed"
# &&
# Note that the Eigen package is installed in a special prefix
# so that we can suppress warnings coming from its include files
# (without suppressing warnings from other include files).
#
# &subhead Debug and Release&&
# If &icode cppad_prefix&& ends in &code cppad_mixed&&,
# all the prefixes must be the same,
# except for the Eigen prefix which has an extra sub-directory at the end.
# In this case,
# &code bin/run_cmake.sh&& will switch between a debug and release
# version of &code cppad_mixed&&
# depending on &cref/build_type/run_cmake.sh/build_type/&&.
# The debug version will have an extra &code .debug&& and the release version
# will have an extra &code .release&& at the end of the prefix.
# In addition, a soft link will be used to map the prefix to either
# the debug or release version.
# Note that the same technique will be used to map the &code build&&
# directory to the debug or release version.
# Also note that if you are using both a debug and release, both versions
# of the  &cref/special requirements/install_unix/Special Requirements/&&
# will need to be installed.
#
# &head cppad_cxx_flags&&
# Extra C++ flags used to configure CppAD
# &codep
cppad_cxx_flags='-std=c++11 -Wall -Wshadow -Wconversion'
# &&
#
# &head cmake_libdir&&
# Sub-directory of each prefix where libraries are installed.
# &codep
cmake_libdir='lib64'
# &&
#
# &head ldlt_cholmod&&
# If yes, use &cref ldlt_cholmod&& LDLT factorization where possible.
# Otherwise always use &cref ldlt_eigen&& for LDLT factorization.
# &codep
ldlt_cholmod='yes'
# &&
#
# &head use_atomic_cholesky&&
# If yes, &code cppad_mixed&& will use
# the &cref sparse_ad_cholesky&& atomic AD operation when computing the
# &cref newton_step&&. Otherwise, an LDLT factorization using
# &code eigen&& , with &code AD<double>&& as the scalar type, is used.
# (Note that the &code cholmod&& LDLT factorization cannot
# be use with and AD scalar type.)
# &codep
use_atomic_cholesky='no'
# &&
#
# &head checkpoint_newton_step&&
# If yes, &code cppad_mixed&& will checkpoint the
# &cref newton_step&&. Otherwise, repeated applications of the Newton step
# are recorded on the AD tape (which should require more memory but may be
# faster).
# &codep
checkpoint_newton_step='no'
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
# &head hide_ipopt_scaling&&
# If yes, the re-scaling done by &code cppad_mixed&& is hidden from ipopt.
# This changes the function and derivative values
# printed during the &cref ipopt_trace&&.  On the other hand, it seems
# that hiding the scaling from Ipopt yields better convergence detection.
# If &icode hide_ipopt_scaling&& is &code no&&,
# the ipopt &code user-scaling&& option is used (ipopt is shown the scaling).
# &codep
hide_ipopt_scaling='yes'
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
# in &icode cppad_cxx_flags&&. Then execute the following commands:
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
	[--use_atomic_cholesky] \\
	[--checkpoint_newton_step] \\
	[--optimize_cppad_function] \\
	[--show_ipopt_scaling] \\
	[--rev_hes_sparsity] \\
	[--release] \\
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
	elif [ "$1" == '--use_atomic_cholesky' ]
	then
		use_atomic_cholesky='yes'
	elif [ "$1" == '--checkpoint_newton_step' ]
	then
		checkpoint_newton_step='yes'
	elif [ "$1" == '--optimize_cppad_function' ]
	then
		optimize_cppad_function='yes'
	elif [ "$1" == '--show_ipopt_scaling' ]
	then
		hide_ipopt_scaling='no'
	elif [ "$1" == '--rev_hes_sparsity' ]
	then
		for_hes_sparsity='no'
	elif [ "$1" == '--release' ]
	then
		build_type='release'
	elif [ "$1" == '--dismod_at_prefix' ]
	then
		cppad_prefix="$HOME/prefix/dismod_at"
		eigen_prefix="$HOME/prefix/dismod_at/eigen"
		ipopt_prefix="$HOME/prefix/dismod_at"
		suitesparse_prefix="$HOME/prefix/dismod_at"
	else
		echo "'$1' is an invalid option"
		bin/run_cmake.sh --help
		exit 1
	fi
	shift
done
# --------------------------------------------------------------------------
if echo "$cppad_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh run_cmake $cppad_prefix $build_type
fi
if echo "$cppad_prefix" | grep '/dismod_at$' > /dev/null
then
	bin/build_type.sh run_cmake $cppad_prefix $build_type
fi
export PKG_CONFIG_PATH=$ipopt_prefix/$cmake_libdir/pkgconfig
# --------------------------------------------------------------------------
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
	-D CMAKE_VERBOSE_MAKEFILE=$verbose_makefile \
	-D CMAKE_BUILD_TYPE=$build_type \
	\
	-D cppad_prefix="$cppad_prefix" \
	-D ipopt_prefix="$ipopt_prefix" \
	-D eigen_prefix="$eigen_prefix" \
	-D suitesparse_prefix="$suitesparse_prefix" \
	\
	-D cppad_cxx_flags="$cppad_cxx_flags" \
	-D cmake_libdir="$cmake_libdir" \
	-D ldlt_cholmod="$ldlt_cholmod" \
	-D use_atomic_cholesky="$use_atomic_cholesky" \
	-D checkpoint_newton_step="$checkpoint_newton_step" \
	-D optimize_cppad_function="$optimize_cppad_function" \
	-D hide_ipopt_scaling="$hide_ipopt_scaling" \
	-D for_hes_sparsity="$for_hes_sparsity" \
	..
# ---------------------------------------------------------------------------
echo 'run_cmake.sh: OK'
exit 0
