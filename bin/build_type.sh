#! /bin/bash -e
# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-22 Bradley M. Bell
# ----------------------------------------------------------------------------
if [ "$0" != "bin/build_type.sh" ] || [ "$3" == '' ]
then
	cat << EOF
bin/build_type.sh program_name install_prefix build_type

program_name:    is the program name for error reporting.
install_prefix:  prefix used to install cppad_mixed
build_type:      this must be either debug or release

1. Sets up following soft link:
	build          -> build.built_type
2. If install_prefix ends in /cppad_mixed, sets up following soft link:
	install_prefix -> install_prefix.build_type
EOF
	exit 1
fi
# -----------------------------------------------------------------------------
kernel=$(uname -s)
if [[ "$kernel" =~ MSYS.* ]]
then
    echo 'Warning: MSYS does not suppor symbolic links'
    exit 0
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
program_name="$1"
install_prefix="$2"
build_type="$3"
if [ "$build_type" != 'debug' ] && [ "$build_type" != 'release' ]
then
	echo "$program_name: build_type is not debug or release"
	exit 1
fi
dollar='$'
total_prefix=`echo $install_prefix | sed \
	-e "s|/cppad_mixed$dollar|/cppad_mixed.$build_type|"`
echo "$total_prefix"
if [ "$total_prefix" == "$install_prefix" ]
then
	echo "$program_name: install prefix does not end in /cppad_mixed"
	echo 'Skipping soft link for install prefix'
else
# -----------------------------------------------------------------------------
	if [ ! -e "$total_prefix" ]
	then
		echo_eval mkdir "$total_prefix"
	fi
	if [ -e "$install_prefix" ]
	then
		if [ ! -L "$install_prefix" ]
		then
			echo "$program_name: $install_prefix is not a symbolic link"
			exit 1
		fi
		echo_eval rm "$install_prefix"
	fi
	echo_eval ln -s $total_prefix $install_prefix
# ---------------------------------------------------------------------------
fi
if [ -e build ]
then
	if [ ! -L build ]
	then
		echo "$program_name: build is not a symbolic link"
		exit 1
	fi
	echo_eval rm build
fi
if [ ! -e build.$build_type ]
then
	echo_eval mkdir build.$build_type
fi
echo_eval ln -s build.$build_type build
# ---------------------------------------------------------------------------
echo 'build_type.sh: OK'
exit 0
