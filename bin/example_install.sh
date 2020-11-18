#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-20 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
#
# $OMhelpKeyCharacter=&
# &begin example_install.sh&& &newlinech #&&
#
# &section An Example Installation&&
#
# &head Syntax&&
# &codei%bin/example_install.sh %existing%&&
#
# &head existing&&
# is either &code replace&& or &code use&&.
# If it is replace, pre-existing installs
# will be replaced (takes more time).
# If it is use, pre-existing installs
# will be used (takes less time).
#
# &head Source&&
# &srcthisfile%0%# BEGIN BASH%# END BASH%1%&&
# &end
#
# BEGIN BASH
if [ $0 != 'bin/example_install.sh' ]
then
	echo 'bin/example_install.sh: must be executed from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# --------------------------------------------------------------------------
if [ "$1" != 'replace' ] && [ "$1" != 'use' ]
then
	echo 'bin/example_install.sh: existing'
	echo 'where existing is replace or use'
	exit 1
fi
existing="$1"
# ---------------------------------------------------------------------------
# set build_type to value in run_cmake.sh
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# set cmake_install_prefix to value in run_cmake.sh
cmd=`grep '^cmake_install_prefix=' bin/run_cmake.sh`
eval $cmd
#
# ipopt_prefix
ipopt_prefix="$cmake_install_prefix"
# --------------------------------------------------------------------------
# remove old version of example_install.log and example_install.err
for ext in log err
do
	if [ -e "example_install.$ext" ]
	then
		echo_eval rm "example_install.$ext"
	fi
done
# --------------------------------------------------------------------------
# set build link to build.debug or build.release depending on build_type
if echo "$cmake_install_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh example_install.sh $cmake_install_prefix $build_type
fi
# -----------------------------------------------------------------------------
# set system_type,
# system_install
# example_install.tmp
if which apt-get >& /dev/null
then
	system_type='debian'
	system_install='sudo apt-get install -y'
	apt list --installed | sed -e 's|  *| |g' > example_install.tmp
elif which dnf >& /dev/null
then
	system_type='red_hat'
	system_install='sudo dnf install -y'
	dnf list installed | sed -e 's|  *| |g' > example_install.tmp
elif which yum >& /dev/null
then
	system_type='red_hat'
	system_install='sudo yum install -y'
	yum list installed | sed -e 's|  *| |g' > example_install.tmp
elif which brew >& /dev/null
then
	system_type='mac_brew'
	system_install='brew install'
	brew list | sed -e 's|  *|\n|g' > example_install.tmp
else
	echo 'Cannot find the system pakcage manager'
	exit 1
fi
# --------------------------------------------------------------------------
# system external installs for normal system requirements
if [ "$system_type" == 'debian' ]
then
	list='
		cmake
		git
		wget
		libblas-dev
		liblapack-dev
		suitesparse-dev
		pkg-config
		g++
		gfortran
		libgsl0-dev
	'
elif [ "$system_type" == 'mac_brew' ]
then
	list='
		cmake
		wget
		lapack
		suite-sparse
		pkg-config
		gcc
		gsl
		openjdk
	'
elif [ "$system_type" == 'red_hat' ]
then
	list='
		cmake
		git
		wget
		blas-devel
		lapack-devel
		suitesparse-devel
		pkgconf
		gcc-c++
		gcc-gfortran
		gsl-devel
		java
	'
	# make sure Ipopt configure sees brew version of javac (not /usr/bin/javac)
	PATH="/usr/local/opt/openjdk/bin:$PATH"
else
	echo 'example_install.sh: script error'
	exit 1
fi
for package in $list
do
	if grep "^$package[^a-zA-Z_]" example_install.tmp > /dev/null
	then
		version=`grep "^$package[^a-zA-Z_]" example_install.tmp | head -1`
		echo "using installed $version"
	elif grep "^$package\$" example_install.tmp > /dev/null
	then
		# brew list case
		echo "using installed $package"
	else
		echo_eval $system_install $package
	fi
done
rm example_install.tmp
# ----------------------------------------------------------------------------
# local external installs for special requirements
for pkg in eigen ipopt cppad
do
	# eval below converts $HOME in $prefix to its value for current user
	case $pkg in
		eigen)
		eval file="$cmake_install_prefix/eigen/include/Eigen/Core"
		;;

		ipopt)
		eval file="$ipopt_prefix/include/coin-or/IpIpoptApplication.hpp"
		;;

		cppad)
		eval file="$cmake_install_prefix/include/cppad/cppad.hpp"
		;;

		*)
		echo 'bin/example_install.sh: program error'
		exit 1
		;;
	esac
	#
	install='true'
	if [ -e "$file" ]
	then
		if [ "$existing" == 'replace' ]
		then
			echo "replacing existing $pkg install"
		else
			install='false'
			echo "using existing $pkg install"
		fi
	fi
	p='example_install'
	if [ "$install" == 'true' ]
	then
		echo "bin/install_$pkg.sh 1>> $p.log 2>> $p.err"
		bin/install_$pkg.sh 1>> $p.log 2>> $p.err
	fi
done
# ----------------------------------------------------------------------------
# cppad_mixed
# ----------------------------------------------------------------------------
# bin/run_cmake.sh
echo "bin/run_cmake.sh 1>> example_install.log 2>> example_install.err"
bin/run_cmake.sh 1>> ../example_install.log 2>> ../example_install.err
#
# change into build directory
echo_eval cd build
#
# make check, speed, and install
for cmd in check speed install
do
	echo "make $cmd 1>> example_install.log 2>> example_install.err"
	make $cmd 1>> ../example_install.log 2>> ../example_install.err
done
cd ..
# ----------------------------------------------------------------------------
echo 'example_install.sh: OK'
exit 0
# END BASH
