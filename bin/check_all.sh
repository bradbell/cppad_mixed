#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/check_all.sh" ]
then
	echo "bin/check_all.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
input=''
while [ "$input" != 'd' ] &&  [ "$input" != 'r' ]
do
	msg='Debug [d], Release [r] ?'
	read -p "$msg" input
done
# -----------------------------------------------------------------------------
list=`ls bin/check_*.sh`
for script in $list
do
	if [ "$script" != 'bin/check_all.sh' ]
	then
		$script
	fi
done
# ----------------------------------------------------------------------------
bin/run_omhelp.sh xml
#
if [ "$input" == 'r' ]
then
	bin/run_cmake.sh --release
else
	bin/run_cmake.sh
fi
#
cd build
make check
make speed
make install
cd ..
bin/check_install.sh
#
echo 'check_all.sh: OK'
