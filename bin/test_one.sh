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
# ---------------------------------------------------------------------------
if [ "$1" == '' ]
then
	echo 'usage: bin/test_one.sh path_to_file'
	exit 1
fi
path_to_file="$1"
dir=`echo $path_to_file | sed -e 's|/.*||'`
name=`echo $path_to_file | sed -e 's|.*/||' -e 's|\.cpp$||'`
if ! grep '^# if 0' $dir/$dir.cpp > /dev/null
then
# ---------------------------------------------------------------------------
cat << EOF > test_one.$$
/This comment expected by bin\\/test_one.sh/b start_run
/This comment also expected by bin\\/test_one.sh/b end_run
b done
#
: start_run
s|^|\\tRUN($name);\\n# if 0\\n|
b done
#
: end_run
s|\$|\\n# endif|
: done
EOF
echo_eval sed -f test_one.$$ -i $dir/$dir.cpp
echo_eval rm test_one.$$
# ---------------------------------------------------------------------------
fi
echo_eval cd build
echo_eval make "check_$dir"
