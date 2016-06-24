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
if [ $0 != 'bin/unix_install.sh' ]
then
	echo 'bin/unix_install.sh: must be executed from its parent directory'
	exit 1
fi
if [ "$1" != 'debug' ] && [ "$1" != 'release' ]
then
	echo 'bin/unix_install.sh: build_type'
	echo 'where build_type is debug or release'
	exit 1
fi
build_type="$1"
# -----------------------------------------------------------------------------
if which apt-get > /dev/null
then
	system_type='debian'
	system_install='sudo apt-get install -y'
elif which dnf > /dev/null
then
	system_type='red_hat'
	system_install='sudo dnf install -y'
elif which yum > /dev/null
then
	system_type='red_hat'
	system_install='sudo yum install -y'
else
	echo 'Cannot find the system pakcage manager'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# --------------------------------------------------------------------------
# system external installs
if [ "$system_type" == 'debian' ]
then
	list='
		cmake
		wget
		libblas-dev
		liblapack-dev
		pkg-config
		g++
		gfortran
		libgsl0-dev
	'
else
	list='
		cmake
		wget
		blas-devel
		lapack-devel
		pkgconfig
		gcc-c++
		gcc-gfortran
		gsl-devel
	'
fi
for package in $list
do
	echo_eval $system_install $package
done
# ----------------------------------------------------------------------------
# local external installs
bin/install_eigen.sh $build_type
bin/install_ipopt.sh $build_type
bin/install_suitesparse.sh $build_type
bin/install_cppad.sh
# ----------------------------------------------------------------------------
# cppad_mixed
# ----------------------------------------------------------------------------
list=`echo $PKG_CONFIG_PATH | sed -e 's|:| |'`
found='no'    # have not yet found ipopt.pc
for dir in $list
do
	if [ -e $dir/ipopt.pc ]
	then
		found='yes'
	fi
done
if [ $found == 'no' ]
then
	echo 'Cannot find ipopt.pc in the PKG_CONFIG_PATH directorys'
	cmd=`grep 'ipopt_prefix='  bin/install_ipopt.sh`
	eval $cmd
	dir=`find $ipopt_prefix -name 'ipopt.pc' | head -1 | sed -e 's|/ipopt.pc||'`
	if [ "$dir" != '' ]
	then
		echo "Execute the following comamnd:"
		if [ "$PKG_CONFIG_PATH" == '' ]
		then
			echo "export PKG_CONFIG_PATH=\"$dir\""
		else
			echo "export PKG_CONFIG_PATH=\"\$PKG_CONFIG_PATH:$dir\""
		fi
	else
		echo 'Perhaps bin/install_ipopt.sh failed ?'
	fi
	exit 1
fi
# ----------------------------------------------------------------------------
#
if [ "$build_type" == 'debug' ]
then
	bin/run_cmake.sh
else
	bin/run_cmake.sh --release
fi
cd build
make check
make speed
make install
cd ..
bin/check_install.sh
# ----------------------------------------------------------------------------
echo 'bin/unix_install.sh: OK'
exit 0
