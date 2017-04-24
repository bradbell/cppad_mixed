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
version='3.2.9'
web_page='https://bitbucket.org/eigen/eigen/get'
# --------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# eigen_prefix
cmd=`grep '^eigen_prefix=' bin/run_cmake.sh`
eval $cmd
# --------------------------------------------------------------------------
cppad_mixed_dir=`echo $eigen_prefix | sed -e 's|/cppad_mixed/.*|/cppad_mixed|'`
if echo "$cppad_mixed_dir" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_eigen $cppad_mixed_dir $build_type
fi
# --------------------------------------------------------------------------
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
cd build/external
# --------------------------------------------------------------------------
if [ ! -e eigen-$version.tar.gz ]
then
	echo_eval wget $web_page/$version.tar.gz
	echo_eval mv $version.tar.gz eigen-$version.tar.gz
fi
if [ -e "eigen-$version" ]
then
	echo_eval rm -rf eigen-$version
fi
#
echo_eval tar -xzf eigen-$version.tar.gz
tar_name=`ls | grep eigen-eigen`
echo_eval mv $tar_name eigen-$version
# --------------------------------------------------------------------------
echo_eval cd eigen-$version
#
echo_eval mkdir build
echo_eval cd build
# --------------------------------------------------------------------------
echo_eval cmake \
	-Wno-dev \
	-DCMAKE_INSTALL_PREFIX=$eigen_prefix \
	-DCMAKE_BUILD_TYPE=$build_type \
	..
echo_eval make install
# --------------------------------------------------------------------------
include_dir="$eigen_prefix/include"
if [ ! -h $include_dir/Eigen ]
then
	echo_eval ln -s $include_dir/eigen3/Eigen $include_dir/Eigen
fi
# -----------------------------------------------------------------------------
echo 'bin/install_eigen.sh: OK'
exit 0
