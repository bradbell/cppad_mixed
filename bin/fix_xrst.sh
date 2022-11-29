#! /bin/bash -e
# $Id:$
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# -----------------------------------------------------------------------------
edit_file()
{  file="$1"
   echo $file
   cp $file ~/trash
   sed -i $file -f temp.sed
   echo '--------------------------------------------------------------------'
   if diff $file ~/trash
   then
      echo "$file: no changes"
   else
      echo "$file: changes above"
   fi
}
# -----------------------------------------------------------------------------
cat << EOF > temp.sed
/check latex in omhelp/! b end
N
N
s/.*/# build developer documentation\\
if [ -e doc ]\n\
then\n\
   echo_eval rm -r doc\n\
fi\n\
echo_eval xrst \\\\\\
   --local_toc \\\\\\
   --html_theme sphinx_rtd_theme \\\\\\
   --output_dir doc/
#
: end
EOF
edit_file bin/check_all.sh
# -----------------------------------------------------------------------------
for name in sphinx/spelling sphinx/preamble.rst
do
   if [ -e $name ]
   then
      echo_eval git rm $name
   fi
done
if [ -e sphinx/rst ]
then
   rm -r sphinx/rst
fi
# -----------------------------------------------------------------------------
cat << EOF > temp.sed
/sphinx\\/spelling/d
s|sphinx/preamble[.]rst|xrst.toml|
s|doc[.]omh|doc.xrst|
s|dev[.]omh|dev.xrst|
EOF
edit_file bin/devel.sh
# -----------------------------------------------------------------------------
echo 'fix_xrst.sh: OK'
exit 0
