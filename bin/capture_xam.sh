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
random_seed='12345'
number_random='35'
quasi_fixed='yes'
trace_optimize_fixed='no'
ipopt_solve='yes'
bool_sparsity='no'
hold_memory='no'
derivative_test='no'
start_near_solution='no'
number_fixed_samples='1000'
number_locations='50'
max_population='25.0'
mean_population='5.0'
mean_logit_probability='-0.25'
std_logit_probability='0.25'
random_constraint='yes'
# ---------------------------------------------------------------------------
if [ "$0" != "bin/capture_xam.sh" ]
then
	echo 'bin/capture_xam.sh: must be executed from its parent directory'
	exit 1
fi
if [ ! -e 'build/speed' ]
then
	echo 'capture_xam.sh: must first run bin/run_cmake.sh'
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
echo_eval make capture_xam
cat << EOF
./capture_xam \
$random_seed \
$number_random \
$quasi_fixed \
$trace_optimize_fixed \
$ipopt_solve \
$bool_sparsity \
$hold_memory \
$derivative_test \
$start_near_solution \
$number_fixed_samples \
$number_locations \
$max_population \
$mean_population \
$mean_logit_probability \
$std_logit_probability \
$random_constraint
EOF
./capture_xam \
	$random_seed \
	$number_random \
	$quasi_fixed \
	$trace_optimize_fixed \
	$ipopt_solve \
	$bool_sparsity \
	$hold_memory \
	$derivative_test \
	$start_near_solution \
	$number_fixed_samples \
	$number_locations \
	$max_population \
	$mean_population \
	$mean_logit_probability \
	$std_logit_probability \
	$random_constraint
# ----------------------------------------------------------------------------
echo 'capture_xam.sh: OK'
exit 0
