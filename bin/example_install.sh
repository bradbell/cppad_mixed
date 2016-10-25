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
# &codei%bin/example_install.sh %build_type%&&
#
# &head build_type&&
# is either $code debug$$ or $code release$$.
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
if [ "$1" != 'debug' ] && [ "$1" != 'release' ]
then
	echo 'bin/example_install.sh: build_type'
	echo 'where build_type is debug or release'
	exit 1
fi
build_type="$1"
# --------------------------------------------------------------------------
for ext in log err
do
	if [ -e "example_install.$ext" ]
	then
		echo_eval rm "example_install.$ext"
	fi
done
# --------------------------------------------------------------------------
prefix="$HOME/prefix/cppad_mixed"
if echo "$prefix" | grep '/cppad_mixed$' > /dev/null
then
	bin/build_type.sh example_install.sh $prefix $build_type
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
	if grep $package example_install.tmp > /dev/null
	then
		version=`grep $package example_install.tmp`
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
		eval file="$prefix/eigen/include/Eigen/Core"
		;;

		ipopt)
		eval file="$prefix/include/coin/IpIpoptApplication.hpp"
		;;

		suitesparse)
		eval file="$prefix/include/cholmod.h"
		;;

		cppad)
		eval file="$prefix/include/cppad/cppad.hpp"
		;;

		*)
		echo 'bin/example_install.sh: program error'
		exit 1
		;;
	esac
	#
	response='none'
	if [ -e "$file" ]
	then
		while [ "$response" != 'r' ] && [ "$response" != 'u' ]
		do
			read -p "replace / use existing $pkg install [r/u] ?" response
		done
	else
		response='r'
	fi
	p='example_install'
	if [ "$response" == 'r' ]
	then
		echo "bin/install_$pkg.sh $build_type 1>> $p.log 2>> $p.err"
		bin/install_$pkg.sh $build_type 1>> $p.log 2>> $p.err
	fi
done
# ----------------------------------------------------------------------------
# cppad_mixed
# ----------------------------------------------------------------------------
list=`echo $PKG_CONFIG_PATH | sed -e 's|:| |g'`
found='no'    # have not yet found ipopt.pc
for dir in $list
do
	echo $dir
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
p='example_install'
cmd='bin/run_cmake.sh'
if [ "$build_type" == 'release' ]
then
	cmd='bin/run_cmake.sh --release'
fi
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
#
sed -i bin/check_install.sh \
	-e "s|^build_type=.*|build_type='$build_type'|"
echo "bin/check_install.sh 1>> $p.log 2>> $p.err"
bin/check_install.sh 1>> $p.log 2>> $p.err
sed -i bin/check_install.sh \
	-e "s|^build_type=.*|build_type='debug'|"
# ----------------------------------------------------------------------------
echo 'bin/example_install.sh: OK'
exit 0
# END BASH
