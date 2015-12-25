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
if [ "$0" != "bin/check_example.sh" ]
then
	echo "bin/check_example.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
dir_list=''
list=`git ls-files example`
for file in $list
do
	name=`echo $file | sed -e 's|example/||'`
	ext=`echo $file | sed -e 's|.*\.||'`
	if [ "$ext" == 'cpp' ]
	then
		if ! grep "$name" example/CMakeLists.txt > /dev/null
		then
			echo "$name is not in example/CMakeLists.txt"
			exit 1
		fi
	fi
	if [ "$ext" == 'cpp' ] && [ "$name" != 'example.cpp' ]
	then
		name=`echo $file | sed -e 's|.*/\([^.]*\)\.cpp$|\1|'`
		if ! grep "RUN($name)" example/example.cpp > /dev/null
		then
			echo "RUN($name) is not in example/example.cpp"
			exit 1
		fi
	fi
done
echo 'check_example.sh: OK'