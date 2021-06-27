#! /bin/bash -e
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_cppad.sh' ]
then
	echo 'bin/install_cppad.sh: must be executed from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# --------------------------------------------------------------------------
# Us same version and hash as in cppad_py.git/bin/get_cppad.sh
web_page='https://github.com/coin-or/CppAD.git'
cppad_version='20210606'
hash_code='6dcc0ad96b239f852f0fb6956eeb9eda806d231b'
# --------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# cmake_install_prefix
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
#
# extra_cxx_flags
cmd=`grep '^extra_cxx_flags=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
#
# cppad_prefix
cppad_prefix="$cmake_install_prefix"
# --------------------------------------------------------------------------
# PKG_CONFIG_PATH
export PKG_CONFIG_PATH=''
for dir in $(find -L $cmake_install_prefix -name 'pkgconfig')
do
    PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$dir"
done
# --------------------------------------------------------------------------
if echo "$cmake_install_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_cppad $cmake_install_prefix $build_type
fi
# --------------------------------------------------------------------------
if [ ! -e external/$build_type ]
then
	mkdir -p external/$build_type
fi
echo_eval cd external/$build_type
# --------------------------------------------------------------------------
if [ ! -e cppad.git ]
then
	echo_eval git clone $web_page cppad.git
fi
pwd
#
echo_eval cd cppad.git
echo_eval git checkout master
echo_eval git pull --ff-only
echo_eval git checkout --quiet $hash_code
check=`grep '^SET(cppad_version' CMakeLists.txt | \
	sed -e 's|^[^"]*"\([^"]*\).*|\1|'`
if [ "$cppad_version" != "$check" ]
then
	echo 'install_cppad.sh: cppad_version number does not agree with hash_code'
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
#
cmake_args="-D CMAKE_VERBOSE_MAKEFILE=0"
cmake_args="$cmake_args -D cppad_prefix=$cppad_prefix"
cmake_args="$cmake_args -D include_ipopt=true"
cmake_args="$cmake_args -D cmake_install_libdirs=$cmake_libdir"
echo "cmake $cmake_args -D cppad_cxx_flags='$extra_cxx_flags' .."
cmake $cmake_args -D cppad_cxx_flags="$extra_cxx_flags" ..
#
echo_eval make install
# -----------------------------------------------------------------------------
echo 'install_cppad.sh: OK'
exit 0
