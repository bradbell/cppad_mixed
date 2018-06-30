#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
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
file_list=`bin/ls_files.sh | sed -n -e '/\.cpp$/p'`
define_list='
	CPPAD_MIXED_VERSION
	CPPAD_MIXED_HAS_SUITESPARSE
	CPPAD_MIXED_LDLT
	CPPAD_MIXED_NULL_PTR
	CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
'
for file in $file_list
do
	# remove comment blocks from file
	sed -e '/\/\*/,/\*\//d' "$file" > temp.$$
	#
	required='no'
	present='no'
	for name in $define_list
	do
		if grep $name temp.$$ > /dev/null
		then
			required='yes'
		fi
	done
	if grep '# *include *<cppad/mixed/configure.hpp>' temp.$$ > /dev/null
	then
		present='yes'
	fi
	# Be over cautious and require a direct include in each file that uses
	# the configure settings.
	if [ "$required" == 'yes' ] && [ "$present" == 'no' ]
	then
		echo "missing: # include <cppad/mixed/configure.hpp>"
		echo "	$file"
		rm temp.$$
		exit 1
	fi
	if [ "$required" == 'no' ] && [ "$present" == 'yes' ]
	then
		echo "unecessary: # include <cppad/mixed/configure.hpp>"
		echo "	$file"
		rm temp.$$
		exit 1
	fi
	rm temp.$$
done
# ---------------------------------------------------------------------------
echo 'check_configure.sh: OK'
exit 0
