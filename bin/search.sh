#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# dismod_at: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$1" == '' ]
then
	echo 'usage: bin/search.sh grep_string'
	exit 1
fi
list=`git ls-files`
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
