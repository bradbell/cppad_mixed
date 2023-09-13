#! /usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-23 Bradley M. Bell
# ----------------------------------------------------------------------------
set -e -u
#
# {xrst_begin example_remove.sh}
# {xrst_comment_ch #}
#
# Undo example_install.sh
# #######################
#
# Syntax
# ******
# | |tab| ``bin/example_remove.sh``
#
# Purpose
# *******
# This will remove any files that were installed by
# :ref:`example_install.sh-name` .
# This assumes that the values of
# :ref:`run_cmake.sh@cmake_install_prefix` and
# :ref:`run_cmake.sh@build_type` have not changed.
#
# Source
# ******
# {xrst_literal
#     # BEGIN BASH
#     # END BASH
# }
#
# {xrst_end example_remove.sh}
# ----------------------------------------------------------------------------
# BEGIN BASH
#
# cmake_install_prefix
eval $(grep '^cmake_install_prefix *=' bin/run_cmake.sh)
#
# build_type
eval $(grep '^build_type *=' bin/run_cmake.sh)
#
# remove cppad_mixed
pushd build
make uninstall
popd
#
# remove externals
external_list=$(find external/$build_type -regex '.*[a-zA-Z_].git')
for repo in $external_list
do
   pushd $repo/build
   make uninstall
   popd
done
#
# remove links
if [ -L $cmake_install_prefix/eigen/include/Eigen ]
then
   rm $cmake_install_prefix/eigen/include/Eigen
fi
#
# prune cmake_install_prefix
echo "Remove empty directories below $cmake_install_prefix"
if [ ! -L $cmake_install_prefix ]
then
   find $cmake_install_prefix -type d -empty -delete
   if [ ! -d $cmake_install_prefix ]
   then
      mkdir $cmake_install_prefix
   fi
else
   find $cmake_install_prefix.$build_type -type d -empty -delete
   if [ ! -d $cmake_install_prefix.$build_type ]
   then
      mkdir $cmake_install_prefix.$build_type
   fi
fi
#
echo 'example_remove.sh: OK'
exit 0
# END BASH
