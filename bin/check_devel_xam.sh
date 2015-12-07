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
if [ "$0" != "bin/check_devel_xam.sh" ]
then
	echo "bin/check_devel_xam.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
dir_list=''
list=`ls example/devel`
for name in $list
do
	if [ -d example/devel/$name ] && [ "$name" != 'new' ]
	then
		dir_list="$dir_list $name"
	fi
done
for dir in $dir_list
do
	list=`ls example/devel/$dir/*.cpp | \
		 sed -e "s|example/devel/$dir/||" -e '/^junk$/d' -e '/^new$/d'`
	for file in $list
	do
		if ! grep "$dir/$file" example/devel/CMakeLists.txt > /dev/null
		then
			echo "$dir/$file is not in example/devel/CMakeLists.txt"
			exit 1
		fi
		name=`echo $file | sed -e 's|\.cpp$||'`
		if ! grep "RUN($name)" example/devel/example_devel.cpp > /dev/null
		then
			echo "RUN($name) is not in example/devel/example_devel.cpp"
			exit 1
		fi
	done
done
echo 'check_devel_xam.sh: OK'
