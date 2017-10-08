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
if [ "$0" != "bin/test_one.sh" ]
then
	echo "bin/test_one.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ---------------------------------------------------------------------------
if [ "$1" == '' ]
then
cat << EOF

usage: bin/test_one.sh option
If option is checkout, a git checkout is run on the files
	example/example.cpp, test_more/test_more.cpp.
Otherwise option is a file name under the example or test_more directory
and the corresponding test is run.

EOF
	exit 1
fi
# ---------------------------------------------------------------------------
if [ "$1" == 'checkout' ]
then
	echo_eval git checkout example/example.cpp
	echo_eval git checkout test_more/test_more.cpp
	echo 'test_one.sh: OK'
	exit 0
fi
file_name="$1"
# ---------------------------------------------------------------------------
found='no'
for find_dir in example test_more
do
	find_path=`find $find_dir -name "$file_name"`
	if [ "$find_path" != '' ]
	then
		test_name=`echo $file_name | sed -e 's|.*/||' -e 's|\.cpp$||'`
		if [ "$find_dir" == 'example' ]
		then
			test_name="${test_name}_xam"
		fi
		found='yes'
cat << EOF > test_one.1
/This comment expected by bin\\/test_one.sh/b start_run
/This comment also expected by bin\\/test_one.sh/b end_run
b done
#
: start_run
s|^|\\tRUN($test_name);\\n# if 0\\n|
b done
#
: end_run
s|\$|\\n# endif|
: done
EOF
		git checkout $find_dir/$find_dir.cpp
		echo_eval sed -f test_one.1 -i $find_dir/$find_dir.cpp
		rm test_one.1
		echo_eval cd build
		echo "make check_$find_dir"
		if ! make "check_$find_dir"
		then
			rm ../test_one.1
			echo "cd build; make check_$find_dir failed"
			echo 'bin/test_one.sh: Error'
			exit 1
		fi
		cd ..
	fi
done
if [ "$found" == 'no' ]
then
	echo "test_one.sh: cannot find $file_name in example or test_more"
	exit 1
fi
# ---------------------------------------------------------------------------
echo 'bin/test_one.sh: OK'
exit 0
