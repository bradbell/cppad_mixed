#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ "$0" != "bin/sort_all.sh" ]
then
   echo "bin/sort_all.sh: must be executed from its parent directory"
   exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# ---------------------------------------------------------------------------
list=$(git grep -l 'BEGIN_SORT_THIS_LINE_PLUS_' | sed -e '/\/sort_all.sh$/d')
for file in $list
do
   sort.sh $file
done
# ---------------------------------------------------------------------------
echo 'sort_all.sh: OK'
exit 0
