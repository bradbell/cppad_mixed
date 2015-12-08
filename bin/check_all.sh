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
if [ "$0" != "bin/check_all.sh" ]
then
	echo "bin/check_all.sh: must be executed from its parent directory"
	exit 1
fi
# ---------------------------------------------------------------------------
bin/check_example.sh
bin/check_include.sh
bin/check_verbatim.sh
bin/check_configure.sh
# ----------------------------------------------------------------------------
bin/run_omhelp.sh xml
#
bin/run_cmake.sh
#
cd build
make check
make speed
make install
cd ..
bin/check_install.sh
#
echo 'check_all.sh: OK'
