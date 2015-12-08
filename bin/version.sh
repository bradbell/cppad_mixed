#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ $0 != 'bin/version.sh' ]
then
	echo 'bin/version.sh: must be executed from its parent directory'
	exit 1
fi
if [ "$1" != 'get' ] && [ "$1" != 'set' ] && [ "$1" != 'copy' ]
then
	echo 'usage: bin/version.sh (get | set | copy) [version]'
	echo 'get:  Gets the current version number from ./CMakeLists.txt'
	echo 'set:  Sets ./CMakeLists.txt version number to version.'
	echo '      If version is not present, uses current yyyymmdd'
	echo 'copy: Copies version number from ./CMakeLists.txt to other files.'
	exit 1
fi
echo_eval() {
     echo $*
     eval $*
}
# -----------------------------------------------------------------------------
if [ "$1" == 'set' ]
then
	if [ "$2" == '' ]
	then
		version=`date +%Y%m%d`
	else
		version="$2"
	fi
	echo $version
	sed -e "s|\(SET(cppad_mixed_version\).*|\1 \"$version\" )|" \
		-i.old  CMakeLists.txt
	if diff CMakeLists.txt.old CMakeLists.txt
	then
		echo 'No change in CMakeLists.txt'
		rm CMakeLists.txt.old
		exit 1
	fi
	rm CMakeLists.txt.old
	echo 'bin/version.sh set: OK'
	exit 0
fi
# -----------------------------------------------------------------------------
# get the current version number
version=`sed -n -e '/SET(cppad_mixed_version/p' CMakeLists.txt | \
	sed -e 's|SET(cppad_mixed_version *"\([0-9]\{8\}\)".*|\1|'`
if ! (echo $version | grep '[0-9]\{8\}') > /dev/null
then
	echo 'bin/version.sh: Cannot find verison number in CMakeLists.txt'
	exit 1
fi
if [ "$1" == 'get' ]
then
	echo "$version"
	exit 0
fi
# -----------------------------------------------------------------------------
# Make the version number in the relevant files is the same
list='
	doc.omh
	omh/install_unix.omh
'
for name in $list
do
	sed -e "s|cppad_mixed-[0-9]\{8\}|cppad_mixed-$version|"  \
		-e "s|version = '[0-9]\{8\}'|version = '$version'|" \
		-i.old $name
	#
	echo '-------------------------------------------------------------'
	echo "diff $name.old $name"
	if diff $name.old $name
	then
		echo '	no difference was found'
	fi
	#
	echo_eval rm $name.old
done
# -----------------------------------------------------------------------------
echo 'bin/version.sh copy: OK'
