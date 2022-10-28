#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
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
