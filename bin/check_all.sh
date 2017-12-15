#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
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
		echo_eval $script
	fi
done
#
# check version number
bin/version.sh check
#
# check latex in omhelp
echo_eval run_omhelp.sh -xml doc
# -----------------------------------------------------------------------------
echo_eval cp bin/run_cmake.sh bin/run_cmake.bak
if [ "$1" == 'r' ]
then
	sed -i bin/run_cmake.sh \
		-e "s|^build_type=.*|build_type='release'|" \
		-e "s|^optimize_cppad_function=.*|optimize_cppad_function='yes'|"
else
	sed -i bin/run_cmake.sh \
		-e "s|^build_type=.*|build_type='debug'|" \
		-e "s|^optimize_cppad_function=.*|optimize_cppad_function='no'|"
fi
if [ "$2" == 'a' ]
then
	sed -i bin/run_cmake.sh \
		-e "s|^use_atomic_cholesky=.*|use_atomic_cholesky='yes'|"
else
	sed -i bin/run_cmake.sh \
		-e "s|^use_atomic_cholesky=.*|use_atomic_cholesky='no'|"
fi
if [ "$2" == 'c' ]
then
	sed -i bin/run_cmake.sh \
		-e "s|^checkpoint_newton_step=.*|checkpoint_newton_step='yes'|"
else
	sed -i bin/run_cmake.sh \
		-e "s|^checkpoint_newton_step=.*|checkpoint_newton_step='no'|"
fi
if diff bin/run_cmake.sh bin/run_cmake.bak > /dev/null
then
	echo_eval rm bin/run_cmake.bak
fi
#
bin/run_cmake.sh
#
cd build
make check
make speed
make install
cd ..
bin/check_install.sh
# -----------------------------------------------------------------------------
if [ -e bin/run_cmake.bak ]
then
	echo_eval mv bin/run_cmake.bak bin/run_cmake.sh
fi
echo 'check_all.sh: OK'
