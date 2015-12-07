#! /usr/bin/python3
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
# 	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
import os
import sys
import sqlite3
# ---------------------------------------------------------------------------
line    = None
connection = None
while line != 'quit' :
	line = input('sqlborwser>')
	start = line.find(':') + 1
	end   = len(line)
	text  = line[start : end ].strip()
	#
	# connect file_name
	if line.startswith('connect:') :
		file_name  = text
		if os.path.isfile(file_name) :
			connection = sqlite3.connect(file_name)
		else :
			print('Error: no such file')
	#
	# sql: command
	elif line.startswith('sql:') :
		command = text
		connection.row_factory = sqlite3.Row
		cursor  = connection.cursor()
		first   = True
		for row in cursor.execute(command) :
			if first :
				print( row.keys() )
				first = False
			print( tuple(row) )
	#
	# quit
	elif line.startswith('quit:') :
		sys.exit(0)
	else :
		print('\tconnect: file_name')
		print('\tsql:     command')
		print('\tquit:')
