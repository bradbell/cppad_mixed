#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/check_whats_new.sh" ]
then
	echo "bin/check_whats_new.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
yy=`date +%y`
file='omh/install_unix.omh'
if ! grep "\$cref whats_new_$yy\$\$ release notes" $file > /dev/null
then
	echo "cannot find: \$cref whats_new_$yy\$\$ release notes"
	echo "in file:     $file"
	exit 1
fi
# -----------------------------------------------------------------------------
echo 'check_whats_new.sh: OK'
exit 0

