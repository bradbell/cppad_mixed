#! /usr/bin/python3
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
#
import subprocess
import sys
import os
# ----------------------------------------------------------------------------
usage = 'bin/capture_xam.py number_runs'
if sys.argv[0] != 'bin/capture_xam.py' or len(sys.argv) != 2 :
	sys.exit(usage)
number_runs = int( sys.argv[1] )
os.chdir('build/speed')
# ----------------------------------------------------------------------------
# captue_xam comand line arugments
random_seed              = '0'
number_fixed_samples     = '1000'
number_locations         = '50'
number_times             = '10'
max_population           = '25'
mean_population          = '5.0'
mean_logit_probability   = '-0.5'
std_logit_probability    = '0.5'
quasi_fixed              = 'true'
random_constraint        = 'true'
trace_optimize_fixed     = 'false'
# ----------------------------------------------------------------------------
# make sure cature_xam is up to date
subprocess.call( [ 'make', 'capture_xam' ] )
# ----------------------------------------------------------------------------
for run_count in range( number_runs ) :
	if run_count % 2 == 0 :
		random_seed       = '0'
		random_constraint = 'false'
	else :
		file_name   = 'capture_xam.' + str(run_count - 1)
		fp          = open(file_name, 'r')
		actual_seed = '0'
		for line in fp :
			if line.startswith('actual_seed') :
				index       = line.find('=') + 1
				actual_seed = line[index:].lstrip()
		assert actual_seed != '0'
		random_seed       = actual_seed
		random_constraint = 'true'
	file_name = 'capture_xam.' + str(run_count)
	fp        = open(file_name, 'w')
	command_list   = [
		'./capture_xam',
		random_seed,
		number_fixed_samples,
		number_locations,
		number_times,
		max_population,
		mean_population,
		mean_logit_probability,
		std_logit_probability,
		quasi_fixed,
		random_constraint,
		trace_optimize_fixed
	]
	separator = ' '
	command   = separator.join(command_list) + ' > ' + file_name
	print(command)
	subprocess.call(command_list, stdout=fp)
