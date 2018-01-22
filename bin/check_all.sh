#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
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
if [ "$2" != 'c' ]  && [ "$2" != 'n' ]
then
	ok='no'
fi
if [ "$ok" != 'yes' ]
then
	echo 'bin/check_all.sh build_type newton_step'
	echo 'build_type  = "d" (debug) or "r" (release)'
	echo 'newton_step = 'c' (checkpointed) or 'n' (no checkpointing)'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
# run bin/check_*.sh
list=`ls bin/check_*.sh`
for script in $list
do
	if [ "$script" != 'bin/check_all.sh' ] \
	&& [ "$script" != 'bin/check_install.sh' ]
	then
		$script
	fi
done
#
# check version number
bin/version.sh check
#
# check latex in omhelp
echo_eval run_omhelp.sh -xml doc
# -----------------------------------------------------------------------------
if ! grep "build_type='debug'" bin/run_cmake.sh > /dev/null
then
	echo 'bin/check_all.sh: bin/run_cmake.sh build_type not debug'
	exit 1
fi
if ! grep "optimize_cppad_function='no'" bin/run_cmake.sh > /dev/null
then
	echo 'bin/check_all.sh: bin/run_cmake.sh optimize_cppad_function not no'
	exit 1
fi
if ! grep "checkpoint_newton_step='no'" bin/run_cmake.sh > /dev/null
then
	echo 'bin/check_all.sh: bin/run_cmake.sh checkpoint_newton_step not no'
	exit 1
fi
if [ "$1" == 'r' ]
then
	sed -i bin/run_cmake.sh \
		-e "s|^build_type=.*|build_type='release'|" \
		-e "s|^optimize_cppad_function=.*|optimize_cppad_function='yes'|"
fi
if [ "$2" == 'c' ]
then
	sed -i bin/run_cmake.sh \
		-e "s|^checkpoint_newton_step=.*|checkpoint_newton_step='yes'|"
fi
# -----------------------------------------------------------------------------
echo 'bin/run_cmake.sh >& cmake.log'
bin/run_cmake.sh >& cmake.log
#
cd build
#
for target in check speed install
do
	echo "make $target >& $target.log"
	make $target >& ../$target.log
done
cd ..
bin/check_install.sh
# -----------------------------------------------------------------------------
for target in cmake check speed install
do
	if grep -i 'warningL:' $target.log
	then
		echo "bin/run_check_all.sh: $target.log is has warnings."
		exit 1
	fi
done
# -----------------------------------------------------------------------------
sed -i bin/run_cmake.sh \
	-e "s|^build_type=.*|build_type='debug'|" \
	-e "s|^optimize_cppad_function=.*|optimize_cppad_function='no'|" \
	-e "s|^checkpoint_newton_step=.*|checkpoint_newton_step='no'|"
# -----------------------------------------------------------------------------
echo 'check_all.sh: OK'
exit 0
