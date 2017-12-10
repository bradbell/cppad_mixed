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
if [ "$0" != 'bin/run_omhelp.sh' ]
then
	echo 'bin/run_omhelp.sh must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
if [ -e 'doc' ]
then
	echo_eval rm -r doc
fi
echo_eval mkdir doc
# -----------------------------------------------------------------------------
cd doc
image_link='https://github.com/bradbell/cppad_mixed'
cmd="omhelp ../doc.omh -debug -noframe -image_link $image_link"
echo "$cmd > omhelp.log"
if ! $cmd  > ../omhelp.log
then
	cat ../omhelp.log
	exit 1
fi
if grep '^OMhelp Warning:' ../omhelp.log
then
	exit 1
fi
# -----------------------------------------------------------------------------
echo 'run_omhelp.sh: OK'
exit 0
