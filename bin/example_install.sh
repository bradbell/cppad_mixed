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
# &srcfile%bin/example_install.sh
#	%0%# BEGIN BASH%# END BASH%1%&&
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
existing="$2"
# ---------------------------------------------------------------------------
# build_type
cmd=`grep '^build_type=' bin/run_cmake.sh`
eval $cmd
#
# cppad_prefix
cmd=`grep '^cppad_prefix=' bin/run_cmake.sh`
eval $cmd
#
# ipopt_prefix
cmd=`grep '^ipopt_prefix=' bin/run_cmake.sh`
eval $cmd
# --------------------------------------------------------------------------
for ext in log err
do
	if [ -e "example_install.$ext" ]
	then
		echo_eval rm "example_install.$ext"
	fi
done
# --------------------------------------------------------------------------
if echo "$cppad_prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh example_install.sh $cppad_prefix $build_type
fi
# -----------------------------------------------------------------------------
if which apt-get >& /dev/null
then
	system_type='debian'
	system_list='apt list --installed'
	system_install='sudo apt-get install -y'
elif which dnf >& /dev/null
then
	system_type='red_hat'
	system_list='dnf list installed'
	system_install='sudo dnf install -y'
elif which yum >& /dev/null
then
	system_type='red_hat'
	system_list='yum list installed'
	system_install='sudo yum install -y'
else
	echo 'Cannot find the system pakcage manager'
	exit 1
fi
# --------------------------------------------------------------------------
# system external installs
$system_list | sed -e 's|  *| |g' > example_install.tmp
if [ "$system_type" == 'debian' ]
then
	list='
		cmake
		git
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
		git
		wget
		blas-devel
		lapack-devel
		pkgconf
		gcc-c++
		gcc-gfortran
		gsl-devel
	'
fi
for package in $list
do
	if grep "^$package[^a-zA-Z_]" example_install.tmp > /dev/null
	then
		version=`grep "^$package[^a-zA-Z_]" example_install.tmp | head -1`
		echo "using $version"
	else
		echo_eval $system_install $package
	fi
done
rm example_install.tmp
# ----------------------------------------------------------------------------
# local external installs
for pkg in eigen ipopt suitesparse cppad
do
	# eval below converts $HOME in $prefix to its value for current user
	case $pkg in
		eigen)
		eval file="$cppad_prefix/eigen/include/Eigen/Core"
		;;

		ipopt)
		eval file="$cppad_prefix/include/coin/IpIpoptApplication.hpp"
		;;

		suitesparse)
		eval file="$cppad_prefix/include/cholmod.h"
		;;

		cppad)
		eval file="$cppad_prefix/include/cppad/cppad.hpp"
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
		echo "bin/install_$pkg.sh $build_type 1>> $p.log 2>> $p.err"
		bin/install_$pkg.sh $build_type 1>> $p.log 2>> $p.err
	fi
done
# ----------------------------------------------------------------------------
# cppad_mixed
# ----------------------------------------------------------------------------
dir=`find -L $ipopt_prefix -name 'ipopt.pc' | sed -e 's|/ipopt.pc||'`
if [ "$dir" == '' ]
then
	echo "Cannot find ipopt.pc in $ipopt_prefix directory"
	exit 1
else
	export PKG_CONFIG_PATH="$dir"
fi
# ----------------------------------------------------------------------------
p='example_install'
cmd='bin/run_cmake.sh'
echo "$cmd 1>> $p.log 2>> $p.err"
$cmd 1>> $p.log 2>> $p.err
#
cd build
for cmd in check speed install
do
	echo "make $cmd 1>> $p.log 2>> $p.err"
	make $cmd 1>> ../$p.log 2>> ../$p.err
done
cd ..
# ----------------------------------------------------------------------------
echo 'bin/example_install.sh: OK'
exit 0
# END BASH
