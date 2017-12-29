#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
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
cd include/cppad/mixed
list=`ls *.hpp`
for file in $list
do
	check=`echo $file | tr 'a-z' 'A-Z'`
	check=`echo $check  | sed -e 's|\.|_|' -e 's|^|CPPAD_MIXED_|'`
	if ! grep "^# ifndef $check" $file > /dev/null
	then
		echo "# ifndef $check"
		echo "missing in file include/cppad/mixed/$file"
		exit 1
	fi
	if ! grep "^# define $check" $file > /dev/null
	then
		echo "# define $check"
		echo "missing in file include/cppad/mixed/$file"
		exit 1
	fi
done
# -----------------------------------------------------------------------------
echo 'check_include.sh: OK'
exit 0

