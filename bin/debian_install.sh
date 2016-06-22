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
if [ $0 != 'bin/debian_install.sh' ]
then
	echo 'bin/debian_install.sh: must be executed from its parent directory'
	exit 1
fi
if ! which apt-get
then
	echo 'bin/debian_install.sh: cannot find apt-get'
	exit 1
fi
if [ "$1" != 'debug' ] && [ "$1" != 'release' ]
then
	echo 'bin/debian_install.sh: build_type'
	echo 'where build_type is debug or release'
	exit 1
fi
build_type="$1"
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# --------------------------------------------------------------------------
# system external installs
list='
	cmake
	libblas-dev
	liblapack-dev
	pkg-config
	g++
	gfortran
	libgsl0-dev
	wget
'
for package in $list
do
	echo_eval sudo apt-get install -y $package
done
# ----------------------------------------------------------------------------
# local external installs
bin/install_eigen.sh $build_type
bin/install_ipopt.sh $build_type
bin/install_suitesparse.sh $build_type
bin/install_cppad.sh
# ----------------------------------------------------------------------------
# cppad_mixed
#
# 2DO: The CMakelists.txt files should be fixed so this is not necessary
list='
	example/CMakeLists.txt
	test_more/CMakeLists.txt
	speed/CMakeLists.txt
	bin/check_install.sh
'
for file in $list
do
	# use .log because it is ignored by git, see ../.gitignore
	backup=`echo $file | sed -e 's|\..*|.log|'`
	if [ ! -e $backup ]
	then
		echo_eval mv $file $backup
		sed -e 's|gslcblas|blas|' $backup > $file
	fi
	if diff $backup $file > /dev/null
	then
		echo "kludge fix for $file failed"
		exit 1
	fi
	echo "$file fixed"
done
#
if [ "$build_type" == 'debug' ]
then
	bin/run_cmakd.sh
else
	bin/run_cmakd.sh --release
fi
cd build
make check
make speed
make install
# ----------------------------------------------------------------------------
echo 'bin/debian_install.sh: OK'
exit 0
