#! /bin/bash -e
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# Parameters described in ar1_xam in documentation.
random_seed='12345'
number_random='10000'
quasi_fixed='yes'
trace_optimize_fixed='no'
ipopt_solve='no'
bool_sparsity='no'
hold_memory='no'
derivative_test='no'
start_near_solution='no'
# ---------------------------------------------------------------------------
if [ "$0" != "bin/ar1_xam.sh" ]
then
	echo 'bin/ar1_xam.sh: must be executed from its parent directory'
	exit 1
fi
if [ ! -e 'build/speed' ]
then
	echo 'ar1_xam.sh: must first run bin/run_cmake.sh'
	exit
fi
#
if [ "$1" != 'normal' ] && [ "$1" != 'callgrind' ] && [ "$1" != 'massif' ]
then
	echo 'usage: bin/ar1_xam.sh (normal|callgrind|massif)'
	exit 1
fi
test_to_run="$1"
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ----------------------------------------------------------------------------
echo_eval cd build/speed
echo_eval make ar1_xam
arguments="
$random_seed
$number_random
$quasi_fixed
$trace_optimize_fixed
$ipopt_solve
$bool_sparsity
$hold_memory
$derivative_test
$start_near_solution
"
arguments=`echo $arguments | sed -e 's|\n| |'`
#
if [ "$test_to_run" == 'normal' ]
then
	echo_eval ./ar1_xam $arguments
fi
if [ "$test_to_run" == 'callgrind' ]
then
	echo_eval valgrind \
		--tool=callgrind \
		--callgrind-out-file=callgrind.out.$$ \
		./ar1_xam $arguments
	echo "view with: kcachegrind build/speed/callgrind.out.$$"
fi
if [ "$test_to_run" == 'massif' ]
then
	echo_eval valgrind \
		--tool=massif \
		--massif-out-file=massif.out.$$ \
		./ar1_xam $arguments
	echo_eval ms_print massif.out.$$ > massif.out
	echo "resutls are in build/speed/massif.out"
fi
# ----------------------------------------------------------------------------
echo 'ar1_xam.sh: OK'
exit 0
