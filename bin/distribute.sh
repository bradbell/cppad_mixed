#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/distribute.sh' ] || [ "$USER" != "bradbell" ]
then
	echo 'This program must be run from its parent directory and'
	echo 'only Brad Bell can run it.'
	exit 1
fi
remote_machine='moby'
remote_directory='/www/moby.ihme.washington.edu/htdocs/bradbell'
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# ---------------------------------------------------------------------------
# create distribution
echo_eval rm -rf doc
echo_eval bin/run_omhelp.sh xml printable
echo_eval bin/run_omhelp.sh xml
echo_eval bin/run_omhelp.sh htm printable
echo_eval bin/run_omhelp.sh htm
#
echo_eval tar -czf doc.tgz doc
# --------------------------------------------------------------------------
# copy to remote machine
echo_eval ssh.sh $remote_machine -p doc.tgz doc.tgz
cat << EOF
------------------------------
Enter the following commands:
	ssh.sh $remote_machine -l
	rm -rf doc
	tar -xzf doc.tgz
	rm -r $remote_directory/cppad_mixed
	cp -r doc $remote_directory/cppad_mixed
	exit
------------------------------
EOF
echo 'distribute.sh: OK'
