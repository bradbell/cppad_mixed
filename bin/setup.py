# $Id:$
#  --------------------------------------------------------------------------
# dismod_at: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------

from distutils.core import setup

setup(
	name = 'dismod_at',
	description = 'Disease Rates as Functions of Age and Time',
	license = 'GNU Affero General Public License version 3.0 or later',
	author  = 'Bradley M. Bell',
	author_email= 'bradbell@uw.edu',
	url='http://moby.ihme.washington.edu/bradbell/dismod_at',
	version = '20151123',
	packages = ['dismod_at' ],
	package_dir = { 'dismod_at' : 'python/dismod_at' },
	scripts     = list(),
	data_files  = list()
)
