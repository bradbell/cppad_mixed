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
def check_value(file, name, value) :
	pattern = '^' + name + " *= *'" + value + "'"
	cmd  = [ 'grep', pattern , 'bin/run_cmake.sh' ]
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
check_value('bin/run_cmake.sh', 'build_type', 'release')
check_value('bin/run_cmake.sh', 'optimize_cppad_function', 'yes')
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
			check_value('bin/run_cmake.sh', 'use_atomic_cholesky',    atomic)
			check_value('bin/run_cmake.sh', 'checkpoint_newton_step', check)
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
