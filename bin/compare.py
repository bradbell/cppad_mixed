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
# $begin compare.py$$ $newlinech #$$
# $spell
#	py
#	cppad
#	cholesky
#	ar1_xam
# $$
# $section Compare Speed for Different Options$$
#
# $head Syntax$$
# $codei%bin/compare.py %program% %use_existing%$$
#
# $head Purpose$$
# This script compares memory size and speed of execution for the following
# cases:
# $icode program$$ equal to
# $cref/ar1_xam/ar1_xam.cpp/$$,
# $cref/capture_xam/capture_xam.cpp/$$,
# $cref/use_atomic_cholesky/run_cmake.sh/use_atomic_cholesky/$$
# equal to $code yes$$, $code no$$, and
# $cref/checkpoint_newton_step/run_cmake.sh/checkpoint_newton_step/$$
# equal to $code yes$$, $code no$$.
#
# $head Program Files$$
# The program files are compiled versions of
# $cref ar1_xam.cpp$$ or $cref capture_xam.cpp$$ and are named
# $codei%
#	build.release/speed/%program%_%atomic%_%checkpoint%
# %$$
# for $icode program$$ equal to $code ar1_xam$$, $code capture_xam$$,
# $icode atomic$$ equal to $code yes$$, $code no$$, and
# $icode checkpoint$$ equal to $code yes$$, $code no$$.
#
# $head Result Files$$
# The result files are the output corresponding to
# $cref ar1_xam.cpp$$ or $cref capture_xam.cpp$$ and are named
# $codei%
#	build.release/speed/%program%_%atomic%_%checkpoint%.out
# %$$
# for $icode program$$ equal to $code ar1_xam$$, $code capture_xam$$,
# $icode atomic$$ equal to $code yes$$, $code no$$, and
# $icode checkpoint$$ equal to $code yes$$, $code no$$.
#
# $head use_existing$$
# This command line argument is either $code yes$$ or $code no$$.
# If it is $code yes$$ then the existing result files are used to
# create the summary report. Otherwise, all of the program files
# and result files are recreated.
#
# $head Summary Report$$
# The following is an example summary report:
# $codei%
# capture_xam                                 (atomic,checkpoint)
#                                (yes,yes)    (yes,no)    (no,yes)     (no,no)
# initialize_kilobytes           124316.73   201595.62    97201.79   201337.34
# initialize_milliseconds          3080.00     4270.00     1620.00     3260.00
# optimize_fixed_milliseconds      7800.00     6420.00     6100.00     4950.00
# information_mat_milliseconds      357.00      449.00      324.00      270.00
# compare.py: OK
# %$$
#
# $head program$$
# For this example case, the command line argument $icode program$$,
# has the value $code capture_xam$$.
#
# $head initialize_kilobytes$$
# The values in this row are
# $cref/initialize_bytes/ar1_xam.cpp/Output/initialize_bytes/$$ divided by 1024
# for the case corresponding to each column.
#
# $head initialize_milliseconds$$
# The values in this row are
# $cref/initialize_seconds/ar1_xam.cpp/Output/initialize_seconds/$$ times 1000
# for the case corresponding to each column.
#
# $head optimize_fixed_milliseconds$$
# The values in this row are
# $cref/optimize_fixed_seconds/ar1_xam.cpp/Output/optimize_fixed_seconds/$$
# times 1000 for the case corresponding to each column.
#
# $head information_mat_milliseconds$$
# The values in this row are
# $cref/information_mat_seconds/ar1_xam.cpp/Output/information_mat_seconds/$$
# times 1000 for the case corresponding to each column.
#
# $head (yes,yes)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = yes%$$ and
# $icode%checkpoint% = yes%$$.
#
# $head (yes,no)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = yes%$$ and
# $icode%checkpoint% = no%$$.
#
# $head (no,yes)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = no%$$ and
# $icode%checkpoint% = yes%$$.
#
# $head (no,no)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = no%$$ and
# $icode%checkpoint% = no%$$.
#
# $end
# ---------------------------------------------------------------------------
import sys
import os
import subprocess
#
usage  = '''usage: bin/compare.py program use_existing
program: is either 'ar1_xam' or 'capture_xam'
use_existing: use existing build.release/speed/*.out files, 'yes' or 'no'
'''
if len(sys.argv) != 3 :
	sys.exit(usage)
if sys.argv[0] != 'bin/compare.py' :
	sys.exit(usage)
#
if sys.argv[1] != 'ar1_xam' and sys.argv[1] != 'capture_xam' :
	sys.exit(usage)
program = sys.argv[1]
#
if sys.argv[2] == 'yes' :
	use_existing = True
elif sys.argv[2] == 'no' :
	use_existing = False
else :
	sys.exit(usage)
# ---------------------------------------------------------------------------

def system_cmd(echo, cmd, file_name) :
	if file_name == None :
		file_out = None
	else :
		file_out = open(file_name, 'w')
	command = ''
	for arg in cmd :
		command += arg + ' '
	if file_out != None :
		command += ' > ' + file_name
	#
	if echo :
		print(command)
	flag = subprocess.call(cmd, stdout = file_out)
	if flag != 0 :
		if not echo :
			print(command)
		sys.exit('compare.py: command failed' )
	if file_out != None :
		file_out.close()
# ---------------------------------------------------------------------------
def check_value(file_name, name, value) :
	pattern = '^' + name + " *= *" + value
	cmd  = [ 'grep', pattern , file_name ]
	echo = False
	system_cmd(echo, cmd, 'compare.tmp')
# ---------------------------------------------------------------------------
# commanmd line options that are the same for both ar1_xam and capture_xam
if program == 'ar1_xam' :
	random_seed='0'
	number_random='1000'
	quasi_fixed='false'
	trace_optimize_fixed='false'
	ipopt_solve='false'
	bool_sparsity='false'
	hold_memory='false'
	derivative_test='false'
	start_near_solution='false'
else :
	random_seed='0'
	number_random='35'
	quasi_fixed='false'
	trace_optimize_fixed='true'
	ipopt_solve='false'
	bool_sparsity='false'
	hold_memory='false'
	derivative_test='false'
	start_near_solution='false'
	number_fixed_samples='1000'
	number_locations='50'
	max_population='25'
	mean_population='5.0'
	mean_logit_probability='-0.5'
	std_logit_probability='0.5'
	random_constraint='true'
# ---------------------------------------------------------------------------
# initialize run_cmake
echo = False
cmd  = [ 'git', 'checkout', 'bin/run_cmake.sh' ]
system_cmd(echo, cmd, None)
cmd  = ['sed' , '-i', 'bin/run_cmake.sh' ]
cmd += [ '-e', 's|^\\(build_type\\)=.*|\\1=\'release\'|' ]
cmd += [ '-e', 's|^\\(optimize_cppad_function\\)=.*|\\1=\'yes\'|' ]
system_cmd(echo, cmd, None)
check_value('bin/run_cmake.sh', 'build_type', "'release'")
check_value('bin/run_cmake.sh', 'optimize_cppad_function', "'yes'")
# ---------------------------------------------------------------------------
for atomic in [ 'no' , 'yes' ] :
	for check in [ 'no' , 'yes' ] :
		program_ac = program + '_' + atomic + '_' + check
		file_name  = 'build.release/speed/' + program_ac + '.out'
		use     = False
		if os.path.isfile( file_name ) :
			use = use_existing
		if not use :
			# set run_cmake.sh for this value of atomic and check
			cmd  = ['sed' , '-i', 'bin/run_cmake.sh' ]
			tmp  = 's|^\\(use_atomic_cholesky\\)=.*|\\1=\'' + atomic + '\'|'
			cmd += [ '-e', tmp ]
			tmp  = 's|^\\(checkpoint_newton_step\\)=.*|\\1=\'' + check + '\'|'
			cmd += [ '-e', tmp ]
			echo = True
			system_cmd(echo, cmd, None)
			#
			value = "'" + atomic + "'"
			check_value('bin/run_cmake.sh', 'use_atomic_cholesky',    value)
			value = "'" + check + "'"
			check_value('bin/run_cmake.sh', 'checkpoint_newton_step', value)
			#
			# execute run_cmake.sh
			cmd = [ 'bin/run_cmake.sh' ]
			echo = True
			system_cmd(echo, cmd, 'run_cmake.log')
			#
			# create build.release/speed/program_ac
			os.chdir('build.release/speed')
			cmd = [ 'make', program ]
			print('make ' + program + ' > make.log')
			echo = False
			system_cmd(echo, cmd, '../../make.log')
			#
			cmd = [ 'mv', program, program_ac ]
			echo = True
			system_cmd(echo, cmd, None)
			#
			# execute program
			cmd = [ './' + program_ac ,
				random_seed ,
				number_random ,
				quasi_fixed   ,
				trace_optimize_fixed ,
				ipopt_solve ,
				bool_sparsity ,
				hold_memory ,
				derivative_test ,
				start_near_solution
			]
			if program == 'capture_xam' :
				cmd += [
					number_fixed_samples,
					number_locations,
					max_population,
					mean_population,
					mean_logit_probability,
					std_logit_probability,
					random_constraint
				]
			file_name =  program_ac + '.out'
			echo      = True
			system_cmd(echo, cmd, file_name)
			#
			os.chdir('../..')
# ---------------------------------------------------------------------------
# restore run_cmake and check_install.sh
echo = False
cmd  = [ 'git', 'checkout', 'bin/run_cmake.sh' ]
system_cmd(echo, cmd, None)
cmd = [ 'git', 'checkout', 'bin/check_install.sh' ]
system_cmd(echo, cmd, None)
# ---------------------------------------------------------------------------
true_or_false = { 'yes':'true', 'no':'false' }
list = [
	('initialize_bytes',         'initialize_kilobytes'),
	('initialize_seconds',       'initialize_milliseconds'),
	('optimize_fixed_seconds',   'optimize_fixed_milliseconds'),
	('information_mat_seconds',  'information_mat_milliseconds')
]
line = '{:28}{:>35s}'.format(program, '(atomic,checkpoint)')
print(line)
line = '{:28s}'.format('')
for atomic in [ 'yes' , 'no' ] :
	for check in [ 'yes' , 'no' ] :
		label = '(' + atomic + ',' + check + ')'
		line += ' {:>11s}'.format(label)
print(line)
for (name,label) in list :
	line = '{:28s}'.format(label)
	for atomic in [ 'yes' , 'no' ] :
		for check in [ 'yes' , 'no' ] :
			file_name  = 'build.release/speed/'
			file_name +=  program + '_' + atomic + '_' + check + '.out'
			#
			value = true_or_false[atomic]
			check_value(file_name, 'use_atomic_cholesky', value)
			value = true_or_false[check]
			check_value(file_name, 'checkpoint_newton_step', value)
			#
			cmd  = [ 'sed', file_name, '-n' ]
			cmd += [ '-e', 's|,||g', '-e', '/^' + name + ' *=/p' ]
			file_name   = 'compare.tmp'
			system_cmd(echo, cmd, file_name)
			#
			file_ptr  = open(file_name, 'r')
			file_data = file_ptr.read()
			exec( file_data )
			value = eval(name)
			if name == 'initialize_bytes' :
				value = value / 1024.0
			else :
				value = value * 1000.0
			line += ' {:11.2f}'.format(value)
	print(line)
cmd = [ 'rm' , 'compare.tmp' ]
system_cmd(echo, cmd, None)
# ---------------------------------------------------------------------------
print('compare.py: OK')
