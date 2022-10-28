# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
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
