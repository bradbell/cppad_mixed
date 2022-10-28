#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ "$0" != "bin/check_if_0.sh" ]
then
   echo "bin/check_if_0.sh: must be executed from its parent directory"
   exit 1
fi
# -----------------------------------------------------------------------------
list=`bin/ls_files.sh`
for file in $list
do
   if grep '^# *if  *0 *$' $file
   then
      echo "is present in $file."
      echo 'Use the following command to fix this ?'
      echo "git checkout $file"
      exit 1
   fi
done
# -----------------------------------------------------------------------------
echo 'check_if_0.sh: OK'
exit 0
