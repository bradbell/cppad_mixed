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
# $begin speed_compare.py$$ $newlinech #$$
# $spell
#	py
#	cppad
#	cholesky
#	ar1_xam
# $$
# $section Compare Speed and Memory for Different Options$$
#
# $head Syntax$$
# $codei%bin/speed_compare.py %option%$$
#
# $head Purpose$$
# This script compares speed of execution and memory size for the following
# cases:
# $codei%
#	%program%      %=% ar1_xam %or% capture_xam
#	%atomic%       %=%      no %or%         yes
#	%checkpoint%   %=%      no %or%         yes
# %$$
# Note that $code ar1_xam$$ refers to $cref ar1_xam.cpp$$,
# $code capture_xam$$ refers to $cref capture_xam.cpp$$,
# $icode atomic$$ refers to
# $cref/use_atomic_cholesky/run_cmake.sh/use_atomic_cholesky/$$,
# and $icode checkpoint$$ refers to
# $cref/checkpoint_newton_step/run_cmake.sh/checkpoint_newton_step/$$.
#
# $head option$$
# This argument is one of the following:
# $code build$$, $code run$$, $code summary$$.
# You must complete a $code build$$ before choosing $code run$$ and
# complete a $code run$$ before choosing $code summary$$.
#
# $head build$$
# If $icode option$$ is $code build$$, the program files are compiled.
# To be specific,
# for $icode program$$ equal to $code ar1_xam$$, $code capture_xam$$,
# $icode atomic$$ equal to $code yes$$, $code no$$, and
# $icode checkpoint$$ equal to $code yes$$, $code no$$,
# $codei%
#	build.release/speed/%program%_%atomic%_%checkpoint%
# %$$
# is created; see $cref ar1_xam.cpp$$ and  $cref capture_xam.cpp$$.
#
# $head run$$
# If $icode option$$ is $code run$$, the output files are created.
# To be specific,
# for $icode program$$ equal to $code ar1_xam$$, $code capture_xam$$,
# $icode atomic$$ equal to $code yes$$, $code no$$, and
# $icode checkpoint$$ equal to $code yes$$, $code no$$,
# $codei%
#	build.release/speed/%program%_%atomic%_%checkpoint%.out
# %$$
# is created; see the heading Output in
# $cref/ar1_xam/ar1_xam.cpp/Output/$$ and
# $cref/capture_xam/capture_xam.cpp/Output/$$.
#
# $head summary$$
# If $icode option$$ is $code summary$$,
# a summary is printed for each program.
# The following is an example summary for $code capture_xam$$:
# $codei%
# capture_xam                                 (atomic, checkpoint)
#                                 (no, no)   (no, yes)   (yes, no)  (yes, yes)
# initialize_kilobytes            98740.49    98740.49    98740.49    98740.49
# initialize_milliseconds           350.00      358.00      340.00      367.00
# optimize_fixed_milliseconds      4100.00     3980.00     4050.00     4110.00
# information_mat_milliseconds     3680.00     1700.00     4220.00     2530.00
# final_kilobytes                281869.35   205188.57   282127.63   205316.74
# %$$
#
# $subhead program$$
# For this summary, $icode program$$ has the value $code capture_xam$$.
#
# $subhead initialize_kilobytes$$
# The values in this row are $icode initialize_bytes$$ divided by 1024
# for the program
# $cref/ar1_xam/ar1_xam.cpp/Output/initialize_bytes/$$ or
# $cref/capture_xam/capture_xam.cpp/Output/initialize_bytes/$$
# for the case corresponding to each column.
#
# $subhead initialize_milliseconds$$
# The values in this row are $icode initialize_seconds$$ times 1000
# for the program
# $cref/ar1_xam/ar1_xam.cpp/Output/initialize_seconds/$$ or
# $cref/capture_xam/capture_xam.cpp/Output/initialize_seconds/$$
# for the case corresponding to each column.
#
# $subhead optimize_fixed_milliseconds$$
# The values in this row are $icode optimize_fixed_seconds$$ times 1000
# for the program
# $cref/ar1_xam/ar1_xam.cpp/Output/optimize_fixed_seconds/$$ or
# $cref/capture_xam/capture_xam.cpp/Output/optimize_fixed_seconds/$$
# for the case corresponding to each column.
#
# $subhead information_mat_milliseconds$$
# The values in this row are $icode information_mat_seconds$$ times 1000
# for the program
# $cref/ar1_xam/ar1_xam.cpp/Output/information_mat_seconds/$$ or
# $cref/capture_xam/capture_xam.cpp/Output/information_mat_seconds/$$
# for the case corresponding to each column.
#
# $subhead final_kilobytes$$
# The values in this row are $icode final_bytes$$ divided by 1024
# for the program
# $cref/ar1_xam/ar1_xam.cpp/Output/final_bytes/$$ or
# $cref/capture_xam/capture_xam.cpp/Output/final_bytes/$$
# for the case corresponding to each column.
#
# $subhead (no, no)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = no%$$ and
# $icode%checkpoint% = no%$$.
#
# $subhead (no, yes)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = no%$$ and
# $icode%checkpoint% = yes%$$.
#
# $subhead (yes, no)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = yes%$$ and
# $icode%checkpoint% = no%$$.
#
# $subhead (yes, yes)$$
# The values in this column correspond to
# $icode%use_atomic_cholesky% = yes%$$ and
# $icode%checkpoint% = yes%$$.
#
# $end
# ---------------------------------------------------------------------------
import sys
import os
import subprocess
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
		sys.exit('speed_compare.py: command failed' )
	if file_out != None :
		file_out.close()
	#
	return
# ---------------------------------------------------------------------------
def check_value(file_name, name, value) :
	pattern = '^' + name + " *= *" + value
	cmd  = [ 'grep', pattern , file_name ]
	echo = False
	system_cmd(echo, cmd, 'compare.tmp')
	os.remove('compare.tmp')
	#
	return
# ---------------------------------------------------------------------------
def quote(value) :
	return "'" + value + "'"
# ---------------------------------------------------------------------------
def build(atomic, checkpoint) :
	# restore run_cmake.sh
	echo = False
	cmd  = [ 'git', 'checkout', 'bin/run_cmake.sh' ]
	system_cmd(echo, cmd, None)
	#
	# modify run_cmake.sh
	cmd  = ['sed' , '-i', 'bin/run_cmake.sh' ]
	cmd += [ '-e', 's|^\\(build_type\\)=.*|\\1=\'release\'|' ]
	cmd += [ '-e', 's|^\\(optimize_cppad_function\\)=.*|\\1=\'yes\'|' ]
	tmp  = 's|^\\(use_atomic_cholesky\\)=.*|\\1=\'' + atomic + '\'|'
	cmd += [ '-e', tmp ]
	tmp  = 's|^\\(checkpoint_newton_step\\)=.*|\\1=\'' + checkpoint + '\'|'
	cmd += [ '-e', tmp ]
	echo = False
	system_cmd(echo, cmd, None)
	#
	# check run_cmake.sh
	check_value('bin/run_cmake.sh', 'build_type', quote('release'))
	check_value('bin/run_cmake.sh', 'optimize_cppad_function', quote('yes'))
	check_value('bin/run_cmake.sh', 'use_atomic_cholesky', quote(atomic))
	check_value(
		'bin/run_cmake.sh', 'checkpoint_newton_step', quote(checkpoint)
	)
	#
	# execute run_cmake.sh
	cmd = [ 'bin/run_cmake.sh' ]
	echo = True
	system_cmd(echo, cmd, 'run_cmake.log')
	#
	# change into the build.release/speed directory
	os.chdir('build.release/speed')
	#
	# build the programs
	for program in [ 'ar1_xam', 'capture_xam' ] :
		program_ac = program + '_' + atomic + '_' + checkpoint
		print('make ' + program_ac + ' > make.log' )
		#
		# build ar1_xam
		cmd  = [ 'make', program ]
		echo = False
		system_cmd(echo, cmd, '../../make.log')
		#
		# move ar1_xam to ar1_xam_atomic_checkpoint
		cmd = [ 'mv', program, program_ac ]
		echo = False
		system_cmd(echo, cmd, None)
	#
	# return to the distribution directory
	os.chdir('../..')
	#
	# restore run_cmake
	cmd  = [ 'git', 'checkout', 'bin/run_cmake.sh' ]
	echo = False
	system_cmd(echo, cmd, None)
	#
	# restore check_install.sh
	cmd = [ 'git', 'checkout', 'bin/check_install.sh' ]
	echo = False
	system_cmd(echo, cmd, None)
	#
	return
# ---------------------------------------------------------------------------
def run(program, atomic, checkpoint) :
	#
	# version of the program that we are running
	program_ac = program + '_' + atomic + '_' + checkpoint
	#
	# change into the build.release/speed directory
	os.chdir('build.release/speed')
	#
	# command line arguments
	if program == 'ar1_xam' :
		random_seed='1234'
		number_random='1000'
		quasi_fixed='yes'
		trace_optimize_fixed='no'
		ipopt_solve='no'
		bool_sparsity='no'
		hold_memory='no'
		derivative_test='no'
		start_near_solution='no'
	else :
		random_seed='5678'
		number_random='35'
		quasi_fixed='yes'
		trace_optimize_fixed='no'
		ipopt_solve='no'
		bool_sparsity='no'
		hold_memory='no'
		derivative_test='no'
		start_near_solution='no'
		number_fixed_samples='1000'
		number_locations='50'
		max_population='25'
		mean_population='5.0'
		mean_logit_probability='-0.25'
		std_logit_probability='0.25'
		random_constraint='yes'
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
	#
	# run program and store output
	file_name =  program_ac + '.out'
	echo      = True
	system_cmd(echo, cmd, file_name)
	#
	# reture to distribution directory
	os.chdir('../..')
	#
	return
# ---------------------------------------------------------------------------
def summary(program) :
	name_list     = [
		('initialize_bytes',         'initialize_kilobytes'),
		('initialize_seconds',       'initialize_milliseconds'),
		('optimize_fixed_seconds',   'optimize_fixed_milliseconds'),
		('information_mat_seconds',  'information_mat_milliseconds'),
		('final_bytes',               'final_kilobytes')
	]
	line = '{:28}{:>36s}'.format(program, '(atomic, checkpoint)')
	print(line)
	line = '{:28s}'.format('')
	for atomic in [ 'no' , 'yes' ] :
		for checkpoint in [ 'no' , 'yes' ] :
			label = '(' + atomic + ', ' + checkpoint + ')'
			line += ' {:>11s}'.format(label)
	print(line)
	for (name,label) in name_list :
		line = '{:28s}'.format(label)
		for atomic in [ 'no' , 'yes' ] :
			for checkpoint in [ 'no' , 'yes' ] :
				program_ac    = program + '_' + atomic + '_' + checkpoint
				file_name  = 'build.release/speed/' + program_ac + '.out'
				#
				check_value(file_name, 'optimize_cppad_function', 'yes')
				check_value(file_name, 'use_atomic_cholesky',      atomic)
				check_value(file_name, 'checkpoint_newton_step',   checkpoint)
				#
				cmd  = [ 'sed', file_name, '-n' ]
				cmd += [ '-e', 's|,||g', '-e', '/^' + name + ' *=/p' ]
				file_name   = 'compare.tmp'
				echo = False
				system_cmd(echo, cmd, file_name)
				#
				file_ptr  = open(file_name, 'r')
				file_data = file_ptr.read()
				exec( file_data )
				value = eval(name)
				if name in [ 'initialize_bytes', 'final_bytes' ] :
					value = value / 1024.0
				else :
					value = value * 1000.0
				line += ' {:11.2f}'.format(value)
		print(line)
	# remove temporary file
	os.remove('compare.tmp')
	return
# ---------------------------------------------------------------------------
# main program
usage  = '''usage: bin/speed_compare.py option
option: is either 'build', 'run', 'summary'
'''
if len(sys.argv) != 2 :
	sys.exit(usage)
if sys.argv[0] != 'bin/speed_compare.py' :
	sys.exit(usage)
#
if not sys.argv[1]  in [ 'build', 'run', 'summary' ] :
	sys.exit(usage)
option = sys.argv[1]
#
if option == 'build' :
	for atomic in [ 'no' , 'yes' ] :
		for checkpoint in [ 'no' , 'yes' ] :
			build(atomic, checkpoint)
elif option == 'run' :
	for program in [ 'ar1_xam', 'capture_xam' ] :
		for atomic in [ 'no' , 'yes' ] :
			for checkpoint in [ 'no' , 'yes' ] :
				run(program, atomic, checkpoint)
else :
	assert option == 'summary'
	for program in [ 'ar1_xam', 'capture_xam' ] :
		summary(program)
		print('')
# ---------------------------------------------------------------------------
print('bin/speed_compare.py ' + sys.argv[1] + ': OK')
