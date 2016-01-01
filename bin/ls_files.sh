#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/ls_files.sh" ]
then
	echo "bin/ls_files.sh: must be executed from its parent directory"
	exit 1
fi
# ---------------------------------------------------------------------------
if [ "$1" != '' ]
then
	echo 'usage: bin/ls_files.sh'
	exit 1
fi
# -----------------------------------------------------------------------------
list_all=`git ls-files`
git ls-files -d "$1" > ls_files.$$
for file in $list_all
do
	pattern=`echo $file | sed -e 's|/|[/]|g' -e 's|^|^|' -e 's|$|$|'`
	if ! grep "$pattern" ls_files.$$ > /dev/null
	then
		echo "$file"
	fi
done
rm ls_files.$$
exit 0
