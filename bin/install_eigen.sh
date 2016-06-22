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
# BEGIN USER_SETTINGS
# Prefix below which eigen will be installed. Use a special sub-directory so
# because warnings are supressed for the directory where eigen is installed.
eigen_prefix="$HOME/prefix/cppad_mixed/eigen"
# END USER_SETTINGS
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_eigen.sh' ]
then
	echo 'bin/install_eigen.sh: must be executed from its parent directory'
	exit 1
fi
if [ "$1" == 'debug' ]
then
	build_type='Debug'
elif [ "$1" == 'release' ]
then
	build_type='Release'
else
	echo 'bin/install_eigen.sh: build_type'
	echo 'where build_type is debug or release'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ---------------------------------------------------------------------------
version='3.2.7'
web_page='https://bitbucket.org/eigen/eigen/get'
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
git_name=`ls | grep eigen-eigen`
echo_eval mv $git_name eigen-$version
# --------------------------------------------------------------------------
echo_eval cd eigen-$version
#
echo 'patch eigen source code'
sed -e 's|Scalar l_ii = 0|Scalar l_ii = Scalar(0)|' \
	-i.old Eigen/src/SparseCore/TriangularSolver.h
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
