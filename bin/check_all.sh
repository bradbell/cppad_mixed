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
if [ "$0" != "bin/check_all.sh" ]
then
	echo 'bin/check_all.sh must be executed from its parent directory'
	exit 1
fi
# ---------------------------------------------------------------------------
if [ "$1" != 'd' ]  &&
   [ "$1" != 'r' ]  &&
   [ "$1" != 'da' ] &&
   [ "$1" != 'ra' ]
then
	echo 'bin/check_all.sh [d|r|da|ra]'
	echo 'd  = debug'
	echo 'da = debug with atomic_cholesky'
	echo 'r  = release'
	echo 'ra = release with atomic_cholesky'
	exit 1
fi
option="$1"
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
list=`ls bin/check_*.sh`
for script in $list
do
	if [ "$script" != 'bin/check_all.sh' ] \
	&& [ "$script" != 'bin/check_install.sh' ]
	then
		$script
	fi
done
# ----------------------------------------------------------------------------
if [ "$option" == 'd' ]
then
	build_type='debug'
	bin/run_cmake.sh
elif [ "$option" == 'da' ]
then
	build_type='debug'
	bin/run_cmake.sh --use_atomic_cholesky
elif [ "$option" == 'r' ]
then
	build_type='release'
	bin/run_cmake.sh --release --optimize_cppad_function
elif [ "$option" == 'ra' ]
then
	build_type='release'
	bin/run_cmake.sh --release --use_atomic_cholesky --optimize_cppad_function
else
	echo 'error in check_all.sh script'
	exit 1
fi
# ----------------------------------------------------------------------------
bin/run_omhelp.sh xml
#
sed -i bin/check_install.sh \
	-e "s|^build_type=.*|build_type='$build_type'|"
#
cd build
make check
make speed
make install
cd ..
bin/check_install.sh
#
if [ "$build_type" == 'release' ]
then
	sed -i bin/check_install.sh \
		-e "s|^build_type=.*|build_type='debug'|"
fi
#
echo 'check_all.sh: OK'
