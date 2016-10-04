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
example/example.cpp and test_more/test_more.cpp.
Otherwise option is a file name within the example or test_more directory
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
# ---------------------------------------------------------------------------
dir=''
for d in example test_more
do
	if [ -e "$d/$1" ]
	then
		dir="$d"
	fi
done
if [ "$dir" == '' ]
then
	echo "test_one.sh: cannot find example/$1 or test_more/$1"
	exit 1
fi
#
file_name="$1"
test_name=`echo $file_name | sed -e 's|.*/||' -e 's|\.cpp$|_xam|'`
# ---------------------------------------------------------------------------
cat << EOF > test_one.$$
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
echo_eval git checkout $dir/$dir.cpp
echo_eval sed -f test_one.$$ -i $dir/$dir.cpp
echo_eval rm test_one.$$
# ---------------------------------------------------------------------------
echo_eval cd build
echo_eval make "check_$dir"
