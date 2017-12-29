#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_ipopt.sh' ]
then
	echo 'bin/install_ipopt.sh: must be executed from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# --------------------------------------------------------------------------
version="Ipopt-3.12.6"
third_party="Mumps Metis"
web_page="http://www.coin-or.org/download/source/Ipopt"
# --------------------------------------------------------------------------
# Get user configuration options from run_cmake.sh
#
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# ipopt_prefix
cmd=`grep '^ipopt_prefix=' bin/run_cmake.sh`
eval $cmd
#
# cmake_libdir
cmd=`grep '^cmake_libdir=' bin/run_cmake.sh`
eval $cmd
# --------------------------------------------------------------------------
# set links for this build_type
if echo "$ipopt_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_ipopt $ipopt_prefix $build_type
fi
export PKG_CONFIG_PATH=$ipopt_prefix/$cmake_libdir/pkgconfig
# --------------------------------------------------------------------------
# change into external directory
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
cd build/external
# --------------------------------------------------------------------------
# get this version of ipopt
if [ ! -e $version.tgz ]
then
	echo_eval wget $web_page/$version.tgz
fi
if [ -e $version ]
then
	echo_eval rm -rf $version
fi
echo_eval tar -xzf $version.tgz
# --------------------------------------------------------------------------
# change into ipopt directory and get third party software
echo_eval cd $version
if [ -e ThirdParty/HSL ]
then
	echo_eval rm -rf ThirdParty/HSL
fi
#
for package in $third_party
do
	echo_eval cd ThirdParty/$package
	echo_eval ./get.$package
	echo_eval cd ../..
done
# ----------------------------------------------------------------------------
# build ipopt
if [ ! -e build ]
then
	echo_eval mkdir build
fi
cd build
# 2DO: remove coin_skip_warn_cxxflags=yes when bug in gcc is fixed; see
# https://github.com/JuliaOpt/Ipopt.jl/issues/13
cat << EOF > config.sh
../configure \\
	$build_type_flag \\
	--enable-shared \\
	--enable-static \\
	--prefix=$ipopt_prefix \\
	--libdir=$ipopt_prefix/$cmake_libdir \\
	--with-blas-lib="-lblas" \\
	--with-lapack-lib="-llapack" \\
	coin_skip_warn_cxxflags=yes
EOF
echo_eval cat config.sh
echo_eval sh config.sh
echo_eval make install
# -----------------------------------------------------------------------------
echo 'install_ipopt.sh: OK'
exit 0
