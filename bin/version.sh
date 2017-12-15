#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ $0 != 'bin/version.sh' ]
then
	echo 'bin/version.sh: must be executed from its parent directory'
	exit 1
fi
file_list='
	doc.omh
'
if [ "$1" == 'get' ] && [ "$2" == '' ]
then
	version.sh get
	exit 0
fi
if [ "$1" == 'copy' ] || [ "$1" == 'check' ]
then
	echo_eval version.sh $1 $file_list
else
	echo_eval version.sh $*
fi
# -----------------------------------------------------------------------------
echo "bin/version.sh $1: OK"
