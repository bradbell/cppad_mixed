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

from distutils.core import setup

setup(
	name = 'cppad_mixed',
	description = 'Disease Rates as Functions of Age and Time',
	license = 'GNU Affero General Public License version 3.0 or later',
	author  = 'Bradley M. Bell',
	author_email= 'bradbell@uw.edu',
	url='http://moby.ihme.washington.edu/bradbell/cppad_mixed',
	version = '20151123',
	packages = ['cppad_mixed' ],
	package_dir = { 'cppad_mixed' : 'python/cppad_mixed' },
	scripts     = list(),
	data_files  = list()
)
