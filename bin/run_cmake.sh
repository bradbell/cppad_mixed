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
#	libdir cholesky
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
# Note that if &icode cppad_prefix&& ends in &code cppad_mixed&&,
# &code run_cmake.sh&& will use a link from the prefix to
# &icode%cppad_prefix%.debug%&& or
# &icode%cppad_prefix%.release%&&
# depending on the choice &icode build_type&&.
# Also note that you can override the default setting (&code debug&&)
# using the command
# &codep
#	bin/run_cmake.sh --release
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
# &head ldlt_cholmod&&
# If YES, use &cref ldlt_cholmod&& LDLT factorization where possible.
# Otherwise always use &cref ldlt_eigen&& for LDLT factorization.
# &codep
ldlt_cholmod='YES'
# &&
#
# &head log_fatal_error&&
# If YES, &code cppad_mixed&& will use
# &cref/fatal_error/public/User Defined Functions/fatal_error/&&
# to report its fatal error messages.
# If NO, fatal errors will be converted to asserts
# (which is useful when running a program in a debugger).
# In addition, &cref/warnings/public/User Defined Functions/warning/&&
# where the context in the debugger is helpful, are also converted to asserts.
# &codep
log_fatal_error='YES'
# &&
#
# &head use_atomic_cholesky&&
# If YES, &code cppad_mixed&& will use
# the &cref sparse_ad_cholesky&& atomic AD operation when computing the
# &cref newton_step&&. Otherwise, an LDLT factorization using
# &code eigen&& , with &code AD<double>&& as the scalar type, is used.
# (Note that the &code cholmod&& LDLT factorization cannot
# be use with and AD scalar type.)
# &codep
use_atomic_cholesky='NO'
# &&
#
# &head checkpoint_newton_step&&
# If YES, &code cppad_mixed&& will checkpoint the
# &cref newton_step&&. Otherwise, repeated applications of the Newton step
# are recorded on the AD tape (which should require more memory but may be
# faster).
# &codep
checkpoint_newton_step='NO'
# &&
#
# &head optimize_cppad_function&&
# If YES, the operation sequence for certain CppAD functions
# will be optimized. This makes the code run faster but in some cases
# it can make debugging more complicated.
# &codep
optimize_cppad_function='NO'
# &&
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
	[--no_log]   \\
	[--use_atomic_cholesky] \\
	[--checkpoint_newton_step] \\
	[--optimize_cppad_function] \\
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
	elif [ "$1" == '--no_log' ]
	then
		log_fatal_error='NO'
	elif [ "$1" == '--use_atomic_cholesky' ]
	then
		use_atomic_cholesky='YES'
	elif [ "$1" == '--checkpoint_newton_step' ]
	then
		checkpoint_newton_step='YES'
	elif [ "$1" == '--optimize_cppad_function' ]
	then
		optimize_cppad_function='YES'
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
# --------------------------------------------------------------------------
if echo "$cppad_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh run_cmake $cppad_prefix $build_type
fi
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
	-D extra_cxx_flags="$extra_cxx_flags" \
	-D cmake_libdir="$cmake_libdir" \
	-D ldlt_cholmod="$ldlt_cholmod" \
	-D log_fatal_error="$log_fatal_error" \
	-D use_atomic_cholesky="$use_atomic_cholesky" \
	-D checkpoint_newton_step="$checkpoint_newton_step" \
	-D optimize_cppad_function="$optimize_cppad_function" \
	..
# ---------------------------------------------------------------------------
echo 'run_cmake.sh: OK'
