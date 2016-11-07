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
# $codei%bin/compare.py %use_existing%$$
#
# $head Purpose$$
# This program compares program memory and speed for different choices
# of two of the $code cppad_mixed$$ configuration options.
# To be specific the
# $cref/use_atomic_cholesky/run_cmake.sh/use_atomic_cholesky/$$ and
# $cref/checkpoint_newton_step/run_cmake.sh/checkpoint_newton_step/$$ options.
#
# $head Program Files$$
# The program files are compiled versions of $cref ar1_xam.cpp$$
# and are named
# $codei%
#	build.release/speed/ar1_xam_%atomic%_%checkpoint%
# %$$
# for $icode atomic$$ equal to $code yes$$, $code no$$ and
# for $icode checkpoint$$ equal to $code yes$$, $code no$$.
# The value of $icode atomic$$ is
# $cref/use_atomic_cholesky/run_cmake.sh/use_atomic_cholesky/$$ and
# $icode checkpoint$$ is
# $cref/checkpoint_newton_step/run_cmake.sh/checkpoint_newton_step/$$ for
# the corresponding results.
#
# $head Result Files$$
# The results files are versions $cref ar1_xam.cpp$$ output and named
# $codei%
#	build.release/speed/ar1_xam_%atomic%_%checkpoint%.out
# %$$
# for $icode atomic$$ equal to $code yes$$, $code no$$ and
# for $icode checkpoint$$ equal to $code yes$$, $code no$$.
# The value of $icode atomic$$ is
# $cref/use_atomic_cholesky/run_cmake.sh/use_atomic_cholesky/$$ and
# $icode checkpoint$$ is
# $cref/checkpoint_newton_step/run_cmake.sh/checkpoint_newton_step/$$ for
# the corresponding results.
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
# (atomic,checkpoint)         (yes,yes)    (yes,no)    (no,yes)     (no,no)
# initialize_kilobytes       4584.84277  5726.60156  3794.27832  4826.60352
# initialize_seconds           86.32590    38.13230     0.17443     0.08411
# optimize_fixed_seconds      231.10520   202.26070     0.28559     0.16530
# information_mat_seconds      29.74650    25.82490     0.00221     0.00170
# %$$
#
# $head initialize_kilobytes$$
# The values in this row are
# $cref/initialize_bytes/ar1_xam.cpp/Output/initialize_bytes/$$ divided by 1024
# for the case corresponding to each column.
#
# $head initialize_seconds$$
# The values in this row are
# $cref/initialize_seconds/ar1_xam.cpp/Output/initialize_seconds/$$
# for the case corresponding to each column.
#
# $head optimize_fixed_seconds$$
# The values in this row are
# $cref/optimize_fixed_seconds/ar1_xam.cpp/Output/optimize_fixed_seconds/$$
# for the case corresponding to each column.
#
# $head information_mat_seconds$$
# The values in this row are
# $cref/information_mat_seconds/ar1_xam.cpp/Output/information_mat_seconds/$$
# for the case corresponding to each column.
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
usage  = '''usage: bin/compare.py use_existing
use_existing: use existing build.release/speed/*.out files, 'yes' or 'no'
'''
if len(sys.argv) != 2 :
	sys.exit(usage)
if sys.argv[0] != 'bin/compare.py' :
	sys.exit(usage)
#
if sys.argv[1] == 'yes' :
	use_existing = True
elif sys.argv[1] == 'no' :
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
# ar1_xam commanmd line options
random_seed='0'
number_random='1000'
quasi_fixed='false'
trace_optimize_fixed='false'
ipopt_solve='false'
bool_sparsity='false'
hold_memory='false'
derivative_test='false'
start_near_solution='false'
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
for atomic in [ 'yes' , 'no' ] :
	for check in [ 'yes' , 'no' ] :
		program    = 'ar1_xam_' + atomic + '_' + check
		file_name  = 'build.release/speed/' + program + '.out'
		use     = False
		if os.path.isfile( file_name ) :
			use = use_existing
		if not use :
			echo = True
			# set run_cmake.sh for this value of atomic and check
			cmd  = ['sed' , '-i', 'bin/run_cmake.sh' ]
			tmp  = 's|^\\(use_atomic_cholesky\\)=.*|\\1=\'' + atomic + '\'|'
			cmd += [ '-e', tmp ]
			tmp  = 's|^\\(checkpoint_newton_step\\)=.*|\\1=\'' + check + '\'|'
			cmd += [ '-e', tmp ]
			system_cmd(echo, cmd, None)
			#
			value = "'" + atomic + "'"
			check_value('bin/run_cmake.sh', 'use_atomic_cholesky',    value)
			value = "'" + check + "'"
			check_value('bin/run_cmake.sh', 'checkpoint_newton_step', value)
			#
			# execute run_cmake.sh
			cmd = [ 'bin/run_cmake.sh' ]
			system_cmd(echo, cmd, 'run_cmake.log')
			#
			# create build.release/speed/program
			os.chdir('build.release/speed')
			cmd = [ 'make', 'ar1_xam' ]
			system_cmd(echo, cmd, 'make.log')
			#
			cmd = [ 'mv', 'ar1_xam', program ]
			system_cmd(echo, cmd, None)
			#
			# execute program
			cmd = [ './' + program ,
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
			file_name =  program + '.out'
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
	'initialize_bytes',
	'initialize_seconds',
	'optimize_fixed_seconds',
	'information_mat_seconds'
]
line = '{:25s}'.format('(atomic,checkpoint)')
for atomic in [ 'yes' , 'no' ] :
	for check in [ 'yes' , 'no' ] :
		label = '(' + atomic + ',' + check + ')'
		line += ' {:>11s}'.format(label)
print(line)
for name in list :
	if name == 'initialize_bytes' :
		line = '{:25s}'.format('initialize_kilobytes')
	else :
		line = '{:25s}'.format(name)
	for atomic in [ 'yes' , 'no' ] :
		for check in [ 'yes' , 'no' ] :
			file_name  = 'build.release/speed/'
			file_name +=  'ar1_xam_' + atomic + '_' + check + '.out'
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
			line += ' {:11.5f}'.format(value)
	print(line)
cmd = [ 'rm' , 'compare.tmp' ]
system_cmd(echo, cmd, None)
# ---------------------------------------------------------------------------
print('compare.py: OK')
