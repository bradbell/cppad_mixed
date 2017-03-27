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
random_seed='0'
number_random='10000'
quasi_fixed='no'
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
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ----------------------------------------------------------------------------
echo_eval cd build/speed
echo_eval make ar1_xam
cat << EOF
./ar1_xam \
$random_seed \
$number_random \
$quasi_fixed \
$trace_optimize_fixed \
$ipopt_solve \
$bool_sparsity \
$hold_memory \
$derivative_test \
$start_near_solution
EOF
./ar1_xam \
	$random_seed \
	$number_random \
	$quasi_fixed \
	$trace_optimize_fixed \
	$ipopt_solve \
	$bool_sparsity \
	$hold_memory \
	$derivative_test \
	$start_near_solution
# ----------------------------------------------------------------------------
echo 'ar1_xam.sh: OK'
exit 0
