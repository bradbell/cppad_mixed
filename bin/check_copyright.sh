#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-14 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/check_copyright.sh" ]
then
	echo "bin/check_copyright.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
list=`git ls-files | sed -e '/^\.gitignore/d'`
for file in $list
do
	text='Copyright (C) 2014-.. University of Washington'
	if ! grep "$text" $file > /dev/null
	then
		echo "copyright missing: use following command to correct this"
		echo "	bin/add_copyright.sh $file"
		exit 1
	fi
done
#
list=`git status | sed -n \
        -e '/^[#\t ]*deleted:/p' \
        -e '/^[#\t ]*modified:/p' \
        -e '/^[#\t ]*both modified:/p' \
        -e '/^[#\t ]*renamed:/p' \
        -e '/^[#\t ]*new file:/p' | \
            sed -e 's/^.*: *//' -e 's/ -> /\n/' | \
                sort -u`
ok='yes'
for file in $list
do
	if [ -e "$file" ] && [ "$file" != '.gitignore' ]
	then
		text='Copyright (C) 2014-15 University of Washington'
		if ! grep "$text" $file > /dev/null
		then
			echo $file
			sed -e 's|Copyright (C) 2014-..|Copyright (C) 2014-15|' -i.$$ $file
			if diff $file.$$ $file
			then
				echo 'bin/check_copyright.sh: program error'
				exit 1
			else
				rm  $file.$$
			fi
			ok='no'
		fi
	fi
done
#
if [ "$ok" == 'no' ]
then
	echo 'bin/check_copyright.sh: Error corrected, run again'
	exit 1
fi
echo 'bin/check_copyright.sh: OK'
