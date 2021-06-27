#! /bin/bash -e
# -----------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# -----------------------------------------------------------------------------
if [ "$0" != "bin/add_copyright.sh" ]
then
	echo "bin/add_copyright.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
if [ "$1" == "" ]
then
	echo "bin/add_copyright.sh: file_name"
	exit 1
fi
file_name="$1"
if [ ! -e "$file_name" ]
then
	echo "$file_name does not exist"
	exit 1
fi
ext=`echo $file_name | sed -e 's/.*\.//'`
if [ "$ext" == "" ]
then
	echo "bin/add_copyright.sh: file_name does not have an extension"
	exit 1
fi
if grep 'Copyright *(C)' $file_name > /dev/null
then
	echo "bin/add_copyright.sh: $file_name already has a copyright message"
	exit 1
fi
# -----------------------------------------------------------------------------
if [ -e add_copyright.$$ ]
then
	echo "rm add_copyright.$$"
	rm bin/add_copyright.$$
fi
# -----------------------------------------------------------------------------
if head -1 $file_name | grep '#!' > /dev/null
then
	specify_script='yes'
	head -1 $file_name > bin/add_copyright.$$
else
	specify_script='no'
	touch bin/add_copyright.$$
fi
# -----------------------------------------------------------------------------
case $ext in
	cpp|hpp)
	cat << EOF  >> bin/add_copyright.$$
// \$Id:\$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
EOF
	;;

	omh)
	cat << EOF  >> bin/add_copyright.$$
\$Id:\$
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
EOF
	;;

	cmake | py | R |  sed | sh | txt )
	cat << EOF  >> bin/add_copyright.$$
# \$Id:\$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
EOF
	;;

	tex )
	cat << EOF  >> bin/add_copyright.$$
% \$Id:\$
%  --------------------------------------------------------------------------
% cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
%           Copyright (C) 2014-21 University of Washington
%              (Bradley M. Bell bradbell@uw.edu)
%
% This program is distributed under the terms of the
%	     GNU Affero General Public License version 3.0 or later
% see http://www.gnu.org/licenses/agpl.txt
% ---------------------------------------------------------------------------
EOF
	;;

	bib )
	cat << EOF  >> bin/add_copyright.$$
@comment{
\$Id:\$
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell<at>uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
}
EOF
	;;

	md)
	cat << EOF  >> bin/add_copyright.$$
<!-- \$Id:\$ -->
# License
<pre>
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
</pre>
EOF
	;;

	*)
	echo "bin/add_copyright.sh: $file_name extension is one of"
	echo 'bib, cpp, hpp, omh, cmake, py, R, sed, sh, tex, txt'
	exit 1
esac
# -----------------------------------------------------------------------------
if [ "$specify_script" == 'yes' ]
then
	sed -n -e '2,$p' $file_name >> bin/add_copyright.$$
else
	cat $file_name >> bin/add_copyright.$$
fi
if [ -x $file_name ]
then
	chmod +x bin/add_copyright.$$
fi
mv bin/add_copyright.$$ $file_name
# -----------------------------------------------------------------------------
echo 'add_copyright.sh: OK'
exit 0
