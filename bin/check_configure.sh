#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ ! -e "bin/check_configure.sh" ]
then
	echo "bin/check_configure.sh: must be executed from its parent directory"
	exit 1
fi
file_list=`git ls-files | sed -n -e '/\.cpp$/p'`
define_list='
	CPPAD_MIXED_VERSION
	CPPAD_MIXED_HAS_SUITESPARSE
	CPPAD_MIXED_NULL_PTR
	CPPAD_MIXED_SET_SPARSITY
'
for file in $file_list
do
	required='no'
	present='no'
	for name in $define_list
	do
		if grep $name $file > /dev/null
		then
			required='yes'
		fi
	done
	if grep '# *include *<cppad_mixed/configure.hpp>' $file > /dev/null
	then
		present='yes'
	fi
	if [ "$required" == 'yes' ] && [ "$present" == 'no' ]
	then
		echo "missing: # include <cppad_mixed/configure.hpp>"
		echo "	$file"
		exit 1
	fi
	if [ "$required" == 'no' ] && [ "$present" == 'yes' ]
	then
		echo "unecessary: # include <cppad_mixed/configure.hpp>"
		echo "	$file"
		exit 1
	fi
done
exit 0
