#! /bin/bash -e
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-20 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# $OMhelpKeyCharacter=&
# &begin capture_xam.sh&& &newlinech #&&
# &spell
#	capture_xam
#	callgrind
#	valgrind
# &&
#
# &section Example Using capture_xam&&
#
# &head Syntax&&
# &codei%bin/capture_xam.sh %test2run%&&
#
# &head test2run&&
# This argument must be one of the following:
#
# &subhead normal&&
# This test will just run &cref/capture_xam/capture_xam.cpp/&&.
#
# &subhead callgrind&&
# This test will run &code capture_xam&& with
# &code valgrind --tool=callgrind&&. This tool does execution profiling.
#
# &subhead massif&&
# This test will run &code capture_xam&& with
# &code valgrind --tool=massif&&. This tool does memory profiling.
#
# &subhead Source Code&&
# &srcthisfile%0%# BEGIN SH%# END SH%1%&&
#
# &end
# ----------------------------------------------------------------------------
# BEGIN SH
random_seed='1234'
number_random='100'
quasi_fixed='no'
trace_optimize_fixed='no'
ipopt_solve='no'
bool_sparsity='no'
hold_memory='no'
derivative_test='no'
start_near_solution='no'
number_fixed_samples='1000'
number_locations='50'
max_population='20.0'
mean_population='5.0'
mean_logit_probability='-0.50'
std_logit_probability='0.25'
random_constraint='yes'
# ---------------------------------------------------------------------------
program='bin/capture_xam.sh'
if [ "$0" != "$program" ]
then
	echo "$program: must be executed from its parent directory"
	exit 1
fi
speed_dir='build/speed'
if [ ! -e "$speed_dir" ]
then
	echo "$program: must first run:"
	echo '	bin/run_cmake.sh'
	exit 1
fi
#
if [ "$1" != 'normal' ] && [ "$1" != 'callgrind' ] && [ "$1" != 'massif' ]
then
	echo "usage: $program (normal|callgrind|massif)"
	exit 1
fi
test2run="$1"
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ----------------------------------------------------------------------------
echo_eval cd $speed_dir
echo_eval make capture_xam
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
$number_fixed_samples
$number_locations
$max_population
$mean_population
$mean_logit_probability
$std_logit_probability
$random_constraint
"
arguments=`echo $arguments | sed -e 's|\n| |'`
#
if [ "$test2run" == 'normal' ]
then
	echo_eval ./capture_xam $arguments
fi
if [ "$test2run" == 'callgrind' ]
then
	echo_eval valgrind \
		--tool=callgrind \
		--callgrind-out-file=callgrind.out.$$ \
		./capture_xam $arguments
	echo "view with: kcachegrind build/speed/callgrind.out.$$"
fi
if [ "$test2run" == 'massif' ]
then
	echo_eval valgrind \
		--tool=massif \
		--massif-out-file=massif.out.$$ \
		./capture_xam $arguments
	echo_eval ms_print massif.out.$$ > massif.out
	echo "resutls are in build/speed/massif.out"
fi
# ----------------------------------------------------------------------------
echo 'capture_xam.sh: OK'
exit 0
# END SH
