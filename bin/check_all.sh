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
ok='yes'
if [ "$1" != 'd' ]  && [ "$1" != 'r' ]
then
	ok='no'
fi
if [ "$2" != 'a' ]  && [ "$2" != 't' ]
then
	ok='no'
fi
if [ "$3" != 'c' ]  && [ "$3" != 'n' ]
then
	ok='no'
fi
if [ "$ok" != 'yes' ]
then
	echo 'bin/check_all.sh build_type cholesky newton_step'
	echo 'build_type  = "d" (debug) or "r" (release)'
	echo 'cholesky    = "a" (atomic) or "t" (taped)'
	echo 'newton_step = 'c' (checkpointed) or 'n' (no checkpointing)'
	exit 1
fi
if [ "$1" == 'r' ]
then
	build_type='release'
	release='--release'
	optimize_cppad_function='--optimize_cppad_function'
else
	build_type='debug'
	release=''
	optimize_cppad_function=''
fi
if [ "$2" == 'a' ]
then
	use_atomic_cholesky='--use_atomic_cholesky'
else
	use_atomic_cholesky=''
fi
if [ "$3" == 'c' ]
then
	checkpoint_newton_step='--checkpoint_newton_step'
else
	checkpoint_newton_step=''
fi
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
bin/run_cmake.sh \
	$release \
	$use_atomic_cholesky \
	$checkpoint_newton_step \
	$optimize_cppad_function
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
