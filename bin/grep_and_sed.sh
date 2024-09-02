# ---------------------------------------------------------------------------
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
# SPDX-FileContributor: 2003-24 Bradley M. Bell
# ---------------------------------------------------------------------------
#
# source bin/grep_and_sed
# Sets the shell variables grep and sed to the gnu version for this system
#
# sed, grep
if [ "$(uname)" == 'Darwin' ]
then
if which gsed > /dev/null
   then
      sed=gsed
   else
      echo 'darwin_names.sh: cannot find gsed (gnu sed) on MacOS system'
      sed=sed
   fi
   if which ggrep > /dev/null
   then
      grep=ggrep
   else
      echo 'darwin_names.sh: cannot find ggrep (gnu grep) on MacOS system'
      grep=grep
   fi
fi
