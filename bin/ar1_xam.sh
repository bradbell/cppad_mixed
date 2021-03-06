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
# $OMhelpKeyCharacter=&
# &begin ar1_xam.sh&& &newlinech #&&
# &spell
#	ar1_xam
#	callgrind
#	valgrind
# &&
#
# &section Example Using ar1_xam&&
#
# &head Syntax&&
# &codei%bin/ar1_xam.sh %test2run%&&
#
# &head test2run&&
# This argument must be one of the following:
#
# &subhead normal&&
# This test will just run &cref/ar1_xam/ar1_xam.cpp/&&.
#
# &subhead callgrind&&
# This test will run &code ar1_xam&& with
# &code valgrind --tool=callgrind&&. This tool does execution profiling.
#
# &subhead massif&&
# This test will run &code ar1_xam&& with
# &code valgrind --tool=massif&&. This tool does memory profiling.
#
# &subhead Source Code&&
# &srcthisfile%0%# BEGIN SH%# END SH%1%&&
#
# &end
# ----------------------------------------------------------------------------
# BEGIN SH
random_seed='123'
number_random='100000'
quasi_fixed='no'
trace_optimize_fixed='no'
ipopt_solve='no'
bool_sparsity='no'
hold_memory='no'
derivative_test='no'
start_near_solution='no'
# ---------------------------------------------------------------------------
program='bin/ar1_xam.sh'
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
"
arguments=`echo $arguments | sed -e 's|\n| |'`
#
if [ "$test2run" == 'normal' ]
then
	echo_eval ./ar1_xam $arguments
fi
if [ "$test2run" == 'callgrind' ]
then
	echo_eval valgrind \
		--tool=callgrind \
		--callgrind-out-file=callgrind.out.$$ \
		./ar1_xam $arguments
	echo "view with: kcachegrind build/speed/callgrind.out.$$"
fi
if [ "$test2run" == 'massif' ]
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
# END SH
