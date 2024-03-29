#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#
# {xrst_begin capture_xam.sh}
# {xrst_spell
#     callgrind
#     valgrind
# }
# {xrst_comment_ch #}
#
# Example Using capture_xam
# #########################
#
# Syntax
# ******
# ``bin/capture_xam.sh`` *test2run*
#
# test2run
# ********
# This argument must be one of the following:
#
# normal
# ======
# This test will just run :ref:`capture_xam<capture_xam.cpp-name>` .
#
# callgrind
# =========
# This test will run ``capture_xam`` with
# ``valgrind --tool=callgrind`` . This tool does execution profiling.
#
# massif
# ======
# This test will run ``capture_xam`` with
# ``valgrind --tool=massif`` . This tool does memory profiling.
#
# Source Code
# ===========
# {xrst_literal
#     BEGIN SH
#     END SH
# }
#
# {xrst_end capture_xam.sh}
# ----------------------------------------------------------------------------
# BEGIN SH
random_seed='123'
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
   echo '  bin/run_cmake.sh'
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
