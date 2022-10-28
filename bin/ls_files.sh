#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ "$0" != "bin/ls_files.sh" ]
then
   echo "bin/ls_files.sh: must be executed from its parent directory"
   exit 1
fi
# ---------------------------------------------------------------------------
if [ "$1" != '' ]
then
   echo 'usage: bin/ls_files.sh'
   exit 1
fi
# -----------------------------------------------------------------------------
list_all=`git ls-files`
git ls-files -d > ls_files.$$
for file in $list_all
do
   pattern=`echo $file | sed -e 's|/|[/]|g' -e 's|^|^|' -e 's|$|$|'`
   if ! grep "$pattern" ls_files.$$ > /dev/null
   then
      echo "$file"
   fi
done
rm ls_files.$$
exit 0
