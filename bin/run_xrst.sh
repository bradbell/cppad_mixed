#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: Bradley M. Bell <bradbell@seanet.com>
# SPDX-FileContributor: 2003-25 Bradley M. Bell
# ----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
   echo $*
   eval $*
}
# -----------------------------------------------------------------------------
if [ "$0" != 'bin/run_xrst.sh' ]
then
   echo 'bin/run_xrst.sh must be run from its parent directory.'
   exit 1
fi
if [ $# == 1 ]
then
   if [ "$1" == --help ]
   then
cat << EOF
bin/run_xrst.sh flags
possible flags
--external_links           check documentation external links
--suppress_spell_warnings  do not check for documentaiton spelling errors
EOF
      exit 0
   fi
fi
#
# extra_flags
extra_flags=''
while [ $# != 0 ]
do
   case "$1" in

      --external_links)
      extra_flags+=" $1 --continue_with_warnings"
      ;;

      --suppress_spell_warnings)
      extra_flags+=" $1"
      ;;

      *)
      echo "bin/run_xrst.sh: command line argument "$1" is not valid"
      exit 1
      ;;

   esac
   #
   shift
done
#
# build
# this will setup the build subdirectory
eval $(grep '^build_type=' bin/run_cmake.sh)
bin/build_type.sh run_xrst $build_type
#
# index_page_name
index_page_name=$(\
   sed -n -e '/^ *--index_page_name*/p' .readthedocs.yaml | \
   sed -e 's|^ *--index_page_name *||' \
)
#
# xrst
echo_eval xrst \
--local_toc \
--target html \
--html_theme sphinx_rtd_theme \
--index_page_name $index_page_name \
--group_list default dev \
$extra_flags
# -----------------------------------------------------------------------------
echo 'run_xrst.sh: OK'
exit 0
