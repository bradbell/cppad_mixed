#! /bin/bash -e
# -----------------------------------------------------------------------------
# This file was written by Bradley M. Bell and is used by multiple packages.
# It has the same license as the other files in this package.
# The following settings are different for each package:
package='cppad_mixed'
files_with_version_number='doc.omh'
# -----------------------------------------------------------------------------
echo_eval() {
     echo $*
     eval $*
}
# -----------------------------------------------------------------------------
ok='yes'
if [ "$1" != 'get' ] \
&& [ "$1" != 'set' ] \
&& [ "$1" != 'date' ] \
&& [ "$1" != 'copy' ] \
&& [ "$1" != 'check' ]
then
	ok='no'
fi
if [ "$1" == 'set' ]
then
	if [ "$2" == '' ]
	then
		ok='no'
	fi
	if [ "$3" != '' ]
	then
		ok='no'
	fi
else
	if [ "$2" != '' ]
	then
		ok='no'
	fi
fi
if [ "$ok" != 'yes' ]
then
	echo "version.sh $*"
	echo 'usage: version.sh (get|date|copy|check)'
	echo '       version.sh set version'
	exit 1
fi
cmd="$1"
# -----------------------------------------------------------------------------
# determine version number
if [ ! -f CMakeLists.txt ]
then
	echo 'version.sh: cannot find ./CMakeLists.txt'
	exit 1
fi
cat << EOF > version.$$
/^SET *( *${package}_version *"[0-9]\{8\}[0-9.]*" *)/! b skip
s|^SET *( *${package}_version *"\\([0-9]\{8\}[0-9.]*\\)" *)|\1|
p
: skip
EOF
# version number in CMakeLists.txt
version=`sed -n -f version.$$ CMakeLists.txt`
rm version.$$
if ! (echo $version | grep '[0-9]\{8\}') > /dev/null
then
	echo "version.sh: Cannot find ${package}_verison number in CMakeLists.txt"
	exit 1
fi
if [ "$cmd" == 'set' ]
then
	# version number on command line
	version="$2"
fi
if [ "$cmd" == 'date' ]
then
	# version number corresponding to current date
	version=`date +%Y%m%d`
fi
# -----------------------------------------------------------------------------
if [ "$cmd" == 'get' ]
then
	echo "$version"
	exit 0
fi
# -----------------------------------------------------------------------------
# cases where are are setting the version
if [ "$cmd" == 'set' ] || [ "$cmd" == 'date' ]
then
cat << EOF > version.$$
s|^SET *( *${package}_version *"[0-9.]*" *)|SET(${package}_version "$version")|
EOF
	sed  -i.old CMakeLists.txt -f version.$$
	rm version.$$
	if diff CMakeLists.txt.old CMakeLists.txt
	then
		echo 'No change to CMakeLists.txt'
	fi
	rm CMakeLists.txt.old
	#
	echo 'version.sh set: OK'
	exit 0
fi
# -----------------------------------------------------------------------------
for file in $files_with_version_number
do
	sed -e "s|$package-[0-9]\\{8\\}[0-9.]*|${package}-$version|" \
		< $file > $file.copy
	if ! diff $file $file.copy > /dev/null
	then
		echo '-------------------------------------------------------------'
		echo "diff $file"
		if diff $file $file.copy
		then
			echo 'version.sh: program error'
			exit 1
		else
			if [ "$cmd" == 'copy' ]
			then
				mv $file.copy $file
			fi
			if [ "$cmd" == 'check' ]
			then
				ok='no'
				rm $file.copy
			fi
		fi
	else
		rm $file.copy
	fi
done
if [ "$ok" != 'yes' ]
then
	echo 'version.sh check: Found differences.'
	exit 1
fi
# ----------------------------------------------------------------------------
if [ "$cmd" != 'copy' ] && [ "$cmd" != 'check' ]
then
	echo 'version.sh: program error'
	exit 1
fi
echo "version.sh $cmd: OK"
exit 0
