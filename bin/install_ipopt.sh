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
version="3.13.2"
coinbrew='https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew'
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
# --------------------------------------------------------------------------
# change into external directory
if [ ! -e build/external ]
then
	mkdir -p build/external
fi
cd build/external
# -----------------------------------------------------------------------------
if [ ! -e coinbrew ]
then
    echo_eval wget $coinbrew
    echo_eval chmod +x coinbrew
fi
if [ ! -e Ipoot ]
then
    ./coinbrew fetch Ipopt@$version --no-prompt
fi
# -----------------------------------------------------------------------------
# klugde necessary until coin or mumps fixes this problem
cat << EOF > junk.f
      program junk
      print*, "Hello World"
      end
EOF
if gfortran -c -fallow-argument-mismatch junk.f >& /dev/null
then
    echo 'Adding -fallow-argument-mismatch to Mumps fortran compiler flags'
    ADD_FCFLAGS='ADD_FCFLAGS=-fallow-argument-mismatch'
else
    ADD_FCFLAGS=''
fi
# -----------------------------------------------------------------------------
echo_eval ./coinbrew build Ipopt@$version \
	--test \
	--no-prompt \
	--verbosity=3 \
    --prefix=$ipopt_prefix \
	--libdir=$ipopt_prefix/$cmake_libdir \
	$ADD_FCFLAGS
#
echo_eval ./coinbrew install Ipopt@$version \
    --no-prompt
# -----------------------------------------------------------------------------
echo 'install_ipopt.sh: OK'
exit 0
