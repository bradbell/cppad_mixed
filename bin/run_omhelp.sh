#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/run_omhelp.sh' ]
then
	echo 'bin/run_omhelp.sh must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
if [ "$1" != 'htm' ] && [ "$1" != 'xml' ] && [ "$1" != 'all' ]
then
	echo "usage: bin/run_omhelp.sh (htm|xml|all)"
	exit 1
fi
option="$1"
# -----------------------------------------------------------------------------
if [ ! -e 'doc' ]
then
	echo_eval mkdir doc
fi
# -----------------------------------------------------------------------------
version=`bin/version.sh get`
rm -rf doc
echo_eval mkdir -p doc/cppad_mixed-$version
for file in `git ls-files`
do
	sub_dir=`echo $file | sed -e 's|/[^/]*$||'`
	if [ "$sub_dir" != "$file" ]
	then
		if [ ! -e doc/cppad_mixed-$version/$sub_dir ]
		then
			mkdir -p doc/cppad_mixed-$version/$sub_dir
		fi
	fi
	cp $file doc/cppad_mixed-$version/$file
done
echo_eval cd doc
echo_eval tar -czf cppad_mixed-$version.tgz cppad_mixed-$version
# -----------------------------------------------------------------------------
#
if [ "$option" == 'htm' ]
then
	flag_list=':'
elif [ "$option" == 'xml' ]
then
	flag_list='-xml:'
else
	flag_list='-xml: -xml:-printable :-printable :'
fi
for pair in $flag_list
do
	flags=`echo $pair | sed -e 's|:| |'`
	cmd="omhelp ../doc.omh -debug -noframe $flags"
	echo "$cmd > omhelp.log"
	if ! $cmd  > ../omhelp.log
	then
		cat ../omhelp.log
		exit 1
	fi
	if grep '^OMhelp Warning:' ../omhelp.log
	then
		exit 1
	fi
done
echo 'run_omhelp.sh: OK'
exit 0
