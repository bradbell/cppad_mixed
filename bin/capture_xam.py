#! /usr/bin/python3
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
#
import subprocess
import sys
import os
import math
import re
# ----------------------------------------------------------------------------
def get_value(file_name, key) :
	fp = open(file_name, 'r')
	for line in fp :
		if line.startswith(key) :
			pair = line.split('=')
			value = pair[1].strip()
			return value
	fp.close()
	# could not find key in the file
	assert False
# ----------------------------------------------------------------------------
usage  = 'bin/capture_xam.py number_runs'
usage += '\nwhere number_runs is an even positive integer'
if sys.argv[0] != 'bin/capture_xam.py' or len(sys.argv) != 2 :
	sys.exit(usage)
number_runs = int( sys.argv[1] )
assert number_runs % 2 == 0
os.chdir('build/speed')
# ----------------------------------------------------------------------------
# captue_xam comand line arugments
# random_seed            = '0'
number_fixed_samples     = '1000'
number_locations         = '50'
number_times             = '10'
max_population           = '50'
mean_population          = '5.0'
mean_logit_probability   = '-0.5'
std_logit_probability    = '0.5'
quasi_fixed              = 'true'
# random_constraint      = 'true'
trace_optimize_fixed     = 'false' # cannot change this setting
# ----------------------------------------------------------------------------
# make sure cature_xam is up to date
subprocess.call( [ 'make', 'capture_xam' ] )
# ----------------------------------------------------------------------------
# create capture_xam.# files
for run_index in range( number_runs ) :
	if run_index % 2 == 0 :
		random_seed       = '0'
		random_constraint = 'false'
	else :
		file_name         = 'capture_xam.' + str(run_index - 1)
		random_seed       = get_value(file_name, 'actual_seed')
		random_constraint = 'true'
	file_name = 'capture_xam.' + str(run_index)
	if os.path.isfile(file_name) :
		print('use existing ' + file_name)
	else :
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
		fp.close()
# ----------------------------------------------------------------------------
# read results into list of dictionarys
result_list = list()
for run_index in range( number_runs ) :
	file_name = 'capture_xam.' + str(run_index)
	fp        = open(file_name, 'r')
	result    = dict()
	for line in fp :
		pair  = line.split('=')
		key   = pair[0].strip()
		value = pair[1].strip()
		result[key] = value
	result_list.append(result)
	fp.close()
# ----------------------------------------------------------------------------
# check that command line arguments are all equal
# except for random_seed and random_constraint
key_list = [
	'number_fixed_samples',
	'number_locations',
	'number_times',
	'max_population',
	'mean_population',
	'mean_logit_probability',
	'std_logit_probability',
	'quasi_fixed',
	'trace_optimize_fixed'
]
for i in range( number_runs ) :
	for key in key_list :
		assert result_list[0][key] == result_list[i][key]
# ----------------------------------------------------------------------------
# check that sequential even and odd files have same actual_seed values
# and different random_constraint values
n_compare = int( number_runs / 2 )
for i in range( n_compare ) :
	even = 2 * i
	odd  = even + 1
	assert result_list[even]['actual_seed'] == result_list[odd]['actual_seed']
	assert result_list[even]['random_constraint'] == 'false'
	assert result_list[odd]['random_constraint'] == 'true'
# ----------------------------------------------------------------------------
# compute and print averages
average_dict = dict()
#
theta_key = [
	'mean_population', 'mean_logit_probability', 'std_logit_probability'
]
time_even      = 0.0
time_odd       = 0.0
theta_sum_even = [ 0.0, 0.0, 0.0 ]
theta_sum_odd  = [ 0.0, 0.0, 0.0 ]
std_sum_even   = [ 0.0, 0.0, 0.0 ]
std_sum_odd    = [ 0.0, 0.0, 0.0 ]
for i in range( n_compare ) :
	even = 2 * i
	odd  = even + 1
	time_even += float( result_list[even]['optimize_fixed_seconds'] )
	time_odd  += float( result_list[odd]['optimize_fixed_seconds'] )
	for j in range(3) :
		theta_true         = result_list[even][ theta_key[j] ]
		theta_estimate     = result_list[even][ theta_key[j] + '_estimate' ]
		theta_std          = result_list[even][ theta_key[j] + '_std' ]
		diff               = abs(float(theta_estimate) - float(theta_true ))
		theta_sum_even[j] += diff;
		std_sum_even[j]   += float( theta_std )
		#
		theta_true         = result_list[odd][ theta_key[j] ]
		theta_estimate     = result_list[odd][ theta_key[j] + '_estimate' ]
		theta_std          = result_list[odd][ theta_key[j] + '_std' ]
		diff               = abs(float(theta_estimate) - float(theta_true ))
		theta_sum_odd[j]  += diff;
		std_sum_odd[j]    += float( theta_std )
fmt      = '{0:35} false = {1:4.2f} , true = {2:4.2f}'
avg_even = time_even / n_compare
avg_odd  = time_odd  / n_compare
line     = fmt.format('optimize_fixed_seconds_avg', avg_even, avg_odd)
average_dict['optimize_fixed_seconds_avg'] = (avg_even, avg_odd)
print(line)
for j in range(3) :
	avg_even = theta_sum_even[j] / n_compare
	avg_odd  = theta_sum_odd[j]  / n_compare
	std_even = std_sum_even[j] / n_compare
	std_odd  = std_sum_odd[j]  / n_compare
	key      = theta_key[j]
	line     = fmt.format(key + '_error_avg', avg_even, avg_odd)
	average_dict[theta_key[j] + '_error_avg'] = (avg_even, avg_odd )
	print(line)
	key      = key + '_std_avg'
	line     = fmt.format(key, std_even, std_odd)
	print(line)
	average_dict[theta_key[j] + '_std_avg'] = (std_even, std_odd )
# ----------------------------------------------------------------------------
# create tex/capture_xam.new
os.chdir('../../tex')
file_name = 'capture_xam.tex'
fp        = open(file_name, 'r')
data      = fp.read()
fp.close()
#
hspace    = '(\\\\hspace\{[0-9.em]*\})'
single    = '( *[-0-9.,e]+ *| *true *| *false *)'
pair      = '( *[-0-9.,e]+ *, *[-0-9.,e]+ *)'
# input singles
for key in result_list[0].keys() :
	value   = result_list[0][key]
	tmp     = '(' + key.replace('_', r'\\_') + ')'
	pattern = '\n' + tmp + '\n' + hspace + '(.*)\n' + hspace + single + '\n'
	replace = '\n\\1\n\\2\\3\n\\4 ' + value + '\n'
	data    = re.sub(pattern, replace, data)
# output pairs
for key in result_list[0].keys() :
	even    = result_list[0][key]
	odd     = result_list[1][key]
	tmp     = '(' + key.replace('_', r'\\_') + ')'
	pattern = '\n' + tmp + '\n' + hspace + pair + '\n'
	replace = '\n\\1\n\\2 ' + even + ' , ' + odd + '\n'
	data    = re.sub(pattern, replace, data)
# average pairs
for key in average_dict.keys() :
	even    = average_dict[key][0]
	even    = '{0:4.2f}'.format(even)
	odd     = average_dict[key][1]
	odd     = '{0:4.2f}'.format(odd)
	tmp     = '(' + key.replace('_', r'\\_') + ')'
	pattern = '\n' + tmp + '\n' + hspace + pair + '\n'
	replace = '\n\\1\n\\2 ' + even + ' , ' + odd + '\n'
	data    = re.sub(pattern, replace, data)
#
file_name = 'capture_xam.new'
fp        = open(file_name, 'w')
fp.write(data)
fp.close()
sys.exit(0)
