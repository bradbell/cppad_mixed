#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/search.sh" ]
then
	echo "bin/search.sh: must be executed from its parent directory"
	exit 1
fi
# ---------------------------------------------------------------------------
if [ "$1" == '' ]
then
	echo 'usage: bin/search.sh grep_string'
	exit 1
fi
list=`bin/ls_files.sh`
for file in $list
do
	# make git has not staged a delete of this file
	if [ -e "$file" ]
	then
		if grep --ignore-case "$1" $file > /dev/null
		then
			echo $file
		fi
	fi
done
