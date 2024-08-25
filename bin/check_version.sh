#! /usr/bin/env bash
set -e -u
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
# SPDX-FileContributor: 2020-24 Bradley M. Bell
# -----------------------------------------------------------------------------
#
# echo_eval
echo_eval() {
   echo $*
   eval $*
}
#
# sed
source bin/grep_and_sed.sh
# -----------------------------------------------------------------------------
# BEGIN: SECTION THAT DEPENDS ON GIT REPOSITORY
#
# stable_version_file
# This file contains the version number unless the branch is master or main.
# The current date is used when the branch is master or main.
stable_version_file='CMakeLists.txt'
#
# version_files
# These are the files that are check to make sure version number is correct.
version_files='
   CMakeLists.txt
   cppad_mixed.xrst
'
#
# create_temp_sed
# The sed script temp.sed is used to check the version number in files above.
function create_temp_sed {
year=$( echo $version | $sed -e 's|\([0-9]\{4\}\).*|\1|' )
cat << EOF > temp.sed
#
# CMakeLists.txt
s|cppad_mixed_version "[0-9]\\{8\\}"|cppad_mixed_version "$version"|
#
# cppad_mixeddd.xrst
s|cppad_mixed_version-[0-9]\\{8\\}|cppad_mixed_version-$version|
EOF
}
# END: SECTION THAT DEPENDS ON GIT REPOSITORY
# -----------------------------------------------------------------------------
if [ $# != 0 ]
then
   echo 'bin/check_version.sh: does not expect any arguments'
   exit 1
fi
if [ "$0" != 'bin/check_version.sh' ]
then
   echo 'bin/check_version.sh: must be executed from its parent directory'
   exit 1
fi
if [ ! -e './.git' ]
then
   echo 'bin/check_version.sh: cannot find ./.git'
   exit 1
fi
# -----------------------------------------------------------------------------
# check_version
check_version() {
   $sed "$1" -f temp.sed > temp.out
   if ! diff "$1" temp.out > /dev/null
   then
      version_ok='no'
      #
      if [ -x "$1" ]
      then
         mv temp.out "$1"
         chmod +x "$1"
      else
         mv temp.out "$1"
      fi
      echo_eval git diff "$1"
   fi
}
#
# branch
branch=$(git rev-parse --abbrev-ref HEAD)
#
# version
version=$(
   $sed -n -e '/cppad_mixd_version=/p' $stable_version_file | \
      $sed -e 's|.*"\([^"]*\)"|\1|'
)
if [ "$branch" == 'master' ] || [ "$branch" == 'main' ]
then
   version=$(date +%Y%m%d )
fi
if echo $branch | grep '^stable/' > /dev/null
then
   if ! echo $version | grep '[0-9]\{4\}[.]0[.]' > /dev/null
   then
      echo 'check_version.sh: stable version does not begin with yyyy.0.'
      exit 1
   fi
fi
#
# version_ok
version_ok='yes'
#
# check_version
create_temp_sed
for file in $version_files
do
   check_version $file
done
#
# ----------------------------------------------------------------------------
if [ "$version_ok" == 'no' ]
then
   echo 'bin/check_version.sh: version numbers were fixed (see above).'
   echo "Re-execute bin/check_version.sh $version ?"
   exit 1
fi
echo 'check_version.sh OK'
exit 0
