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
if [ "$0" != "bin/check_include.sh" ]
then
	echo "bin/check_include.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
cd include/dismod_at
list=`ls *.hpp`
for file in $list
do
	check=`echo $file | tr 'a-z' 'A-Z'`
	check=`echo $check  | sed -e 's|\.|_|' -e 's|^|DISMOD_AT_|'`
	if ! grep "^# ifndef $check" $file > /dev/null
	then
		echo "# ifndef $check"
		echo "missing in file include/dismod_at/$file"
		exit 1
	fi
	if ! grep "^# define $check" $file > /dev/null
	then
		echo "# define $check"
		echo "missing in file include/dismod_at/$file"
		exit 1
	fi
done
echo 'bin/check_include.sh: OK'

