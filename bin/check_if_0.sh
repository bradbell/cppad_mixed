#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/check_if_0.sh" ]
then
	echo "bin/check_if_0.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
list=`bin/ls_files.sh`
for file in $list
do
	if grep '^# *if  *0 *$' $file
	then
		echo "is present in $file."
		echo 'Use the following command to fix this ?'
		echo "git checkout $file"
		exit 1
	fi
done
echo 'bin/check_if_0.sh: OK'

