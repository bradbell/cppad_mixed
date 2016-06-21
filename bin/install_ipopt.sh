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
# prefix below which ipopt will be installed
ipopt_prefix="$HOME/prefix/cppad_mixed"
# END USER_SETTINGS
# ---------------------------------------------------------------------------
if [ $0 != 'bin/install_ipopt.sh' ]
then
	echo 'bin/install_ipopt.sh: must be executed from its parent directory'
	exit 1
fi
if [ "$1" == 'debug' ]
then
	build_type='--enable-debug'
elif [ "$1" == 'release' ]
then
	build_type=''
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
version="Ipopt-3.11.9"
third_party="Mumps Metis"
web_page="http://www.coin-or.org/download/source/Ipopt"
# --------------------------------------------------------------------------
if [ -e /usr/lib64 ]
then
	libdir='lib64'
else
	libdir='lib'
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
	echo_eval curl -O $web_page/$version.tgz
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
cat << EOF > config.sh
../configure \\
	$build_type \\
	--prefix=$ipopt_prefix \\
	--libdir=$ipopt_prefix/$libdir \\
	--with-blas-lib="-lblas" \\
	--with-lapack-lib="-llapack"
EOF
echo_eval cat config.sh
echo_eval sh config.sh
echo_eval make install | tee make.log
