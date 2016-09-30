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
#
# Prefix below which ipopt will be installed.
# If this directory ends with /cppad_mixed, separate directories are used
# for the debug and release versions.
ipopt_prefix="$HOME/prefix/cppad_mixed"
# END USER_SETTINGS
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_ipopt.sh' ]
then
	echo 'bin/install_ipopt.sh: must be executed from its parent directory'
	exit 1
fi
build_type="$1"
if [ "$build_type" == 'debug' ]
then
	build_type_flag='--enable-debug'
elif [ "$build_type" == 'release' ]
then
	build_type_flag=''
else
	echo 'bin/install_ipopt.sh: build_type'
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
version="Ipopt-3.12.6"
third_party="Mumps Metis"
web_page="http://www.coin-or.org/download/source/Ipopt"
libdir=`bin/libdir.sh`
# --------------------------------------------------------------------------
if echo "$ipopt_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh install_ipopt $ipopt_prefix $build_type
fi
export PKG_CONFIG_PATH=$ipopt_prefix/$libdir/pkgconfig
# --------------------------------------------------------------------------
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
cd build/external
# --------------------------------------------------------------------------
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
	--prefix=$ipopt_prefix \\
	--libdir=$ipopt_prefix/$libdir \\
	--with-blas-lib="-lblas" \\
	--with-lapack-lib="-llapack" \\
	coin_skip_warn_cxxflags=yes
EOF
echo_eval cat config.sh
echo_eval sh config.sh
echo_eval make install | tee make.log
# -----------------------------------------------------------------------------
echo 'bin/install_ipopt.sh: OK'
exit 0
