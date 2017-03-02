#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_cppad.sh' ]
then
	echo 'bin/install_cppad.sh: must be executed from its parent directory'
	exit 1
fi
build_type="$1"
if [ "$build_type" != 'debug' ] && [ "$build_type" != 'release' ]
then
	echo 'bin/install_cppad.sh: build_type'
	echo 'where build_type is debug or release'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# --------------------------------------------------------------------------
web_page='https://github.com/coin-or/CppAD.git'
hash_key='8c8a107cda142a48d92e00d8249b8ced812ea772'
version='20170226'
# --------------------------------------------------------------------------
# cppad_prefix
cmd=`grep '^cppad_prefix=' bin/run_cmake.sh`
eval $cmd
#
# cppad_cxx_flags
cmd=`grep '^cppad_cxx_flags=' bin/run_cmake.sh`
eval $cmd
# --------------------------------------------------------------------------
if echo "$cppad_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_cppad $cppad_prefix $build_type
fi
# --------------------------------------------------------------------------
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
echo_eval cd build/external
# --------------------------------------------------------------------------
if [ ! -e cppad-$version ]
then
	echo_eval git clone $web_page cppad-$version
fi
pwd
#
echo_eval cd cppad-$version
echo_eval git checkout --quiet $hash_key
check=`bin/version.sh get`
if [ "$version" != "$check" ]
then
	echo 'version number does not agree with hash_key'
	exit 1
fi
if ! bin/version.sh check
then
	echo 'version in CMakeLists.txt does not agree with rest of source'
	exit 1
fi
# -----------------------------------------------------------------------------
#
if [ ! -e build ]
then
	mkdir build
fi
echo_eval cd build
if [ -e CMakeCache.txt ]
then
	rm CMakeCache.txt
fi
# -----------------------------------------------------------------------------
if [ -e /usr/lib64 ]
then
	libdirs="'lib64;lib'"
else
	libdirs="'lib;lib64'"
fi
#
cmake_args="-D CMAKE_VERBOSE_MAKEFILE=0"
cmake_args="$cmake_args -D cppad_prefix=$cppad_prefix"
cmake_args="$cmake_args -D cmake_install_libdirs=$libdirs"
echo "cmake $cmake_args -D cppad_cxx_flags='$cppad_cxx_flags' .."
cmake $cmake_args -D cppad_cxx_flags="$cppad_cxx_flags" ..
#
echo_eval make install
# -----------------------------------------------------------------------------
echo 'bin/install_cppad.sh: OK'
exit 0
