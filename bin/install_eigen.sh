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
if [ $0 != 'bin/install_eigen.sh' ]
then
	echo 'bin/install_eigen.sh: must be executed from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ---------------------------------------------------------------------------
version='3.3.7'
web_page='https://gitlab.com/libeigen/eigen.git'
# --------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
#
# cmake_install_prefix
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
#
# eigen_prefix
eigen_prefix="$cmake_install_prefix/eigen"
# --------------------------------------------------------------------------
if echo "$cmake_install_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_eigen $cmake_install_prefix $build_type
fi
# --------------------------------------------------------------------------
if [ ! -e external/$build_type ]
then
	mkdir -p external/$build_type
fi
echo_eval cd external/$build_type
# --------------------------------------------------------------------------
if [ ! -e eigen.git ]
then
	echo_eval git clone $web_page eigen.git
fi
echo_eval cd eigen.git
echo_eval git checkout master
echo_eval git pull --ff-only
echo_eval git checkout $version
#
if [ ! -e build ]
then
	echo_eval mkdir build
fi
echo_eval cd build
# --------------------------------------------------------------------------
echo_eval cmake \
	-Wno-dev \
	-D CMAKE_INSTALL_PREFIX="$eigen_prefix" \
	-D CMAKE_BUILD_TYPE="$build_type" \
	-D PKGCONFIG_INSTALL_DIR="$cmake_install_prefix/$cmake_libdir/pkgconfig" \
	..
echo_eval make install
# --------------------------------------------------------------------------
include_dir="$eigen_prefix/include"
if [ ! -h $include_dir/Eigen ]
then
	echo_eval ln -s $include_dir/eigen3/Eigen $include_dir/Eigen
fi
# -----------------------------------------------------------------------------
echo 'install_eigen.sh: OK'
exit 0
