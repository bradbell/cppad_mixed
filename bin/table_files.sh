#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/table_files.sh' ]
then
	echo 'bin/table_files.sh must be executed from its parent directory'
	exit 1
fi
if [ "$1" == '' ]
then
	echo 'usage: bin/table_files.sh table_name'
	echo 'lists the files that change when a table changes.'
	exit 1
fi
name="$1"
# ---------------------------------------------------------------------------
# 7 is lenght of ${name} in text below
wgt=''
wlt=''
len=`expr length $name`
if [ $len -gt 7 ]
then
	diff=`expr $len - 7`
	for i in $(seq 1 $diff)
	do
		wgt="$wgt "
	done
fi
if [ $len -lt 7 ]
then
	diff=`expr 7 - $len`
	for i in $(seq 1 $diff)
	do
		wlt="$wlt "
	done
fi
cat << EOF
omh/table/${name}_table.omh$wlt                   user doc for $name table
include/cppad_mixed/get_${name}_table.hpp$wlt       C++ struct for $name table
devel/table/get_${name}_table.cpp$wlt             read $name table in C++
python/cppad_mixed/create_database.py$wgt           create all tables
example/table/${name}_table.py$wlt                $name table example
example/table/create_database.py$wgt              create_database.py example
example/get_started/get_started_db.py$wgt         calls create_database
example/user/*.py$wgt                             calls create_database
example/devel/table/get_${name}_table_xam.cpp$wlt C++ example for $name table
test/user/*.py$wgt                                calls create_database
EOF
