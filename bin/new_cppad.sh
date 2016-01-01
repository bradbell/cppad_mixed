# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/new_cppad.sh' ]
then
	echo 'bin/new_cppad.sh must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
dir="$HOME/prefix/cppad_mixed/include/cppad"
echo_eval rm -r $dir
echo_eval cp -r $HOME/cppad.git/cppad $dir
echo "sed -e 's|HAS_COLPACK .*|HAS_COLPACK 0|' -i $dir/configure.hpp"
sed -e 's|HAS_COLPACK .*|HAS_COLPACK 0|' -i $dir/configure.hpp
