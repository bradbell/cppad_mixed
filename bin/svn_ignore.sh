# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-14 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
# 	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# -------------------------------------------------------------------------- */
if [ $0 != "bin/svn_ignore.sh" ]
then
	echo "bin/svn_ignore.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
cat << EOF > bin/svn_ignore.$$
MANIFEST
dist
build
doc
doc.tgz
example.db
junk
junk.*
*.log
*.pyc
new
temp
temp.sh
EOF
svn propset svn:ignore --recursive -F bin/svn_ignore.$$ .
rm bin/svn_ignore.$$
#
exit 0
