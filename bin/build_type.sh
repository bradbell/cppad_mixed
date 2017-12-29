#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/build_type.sh" ] || [ "$3" == '' ]
then
	cat << EOF
bin/build_type.sh program_name install_prefix build_type
sets up soft links for the following directories:
	install_prefix -> install_prefix.build_type
	build          -> build.built_type
program_name:    is the program name for error reporting.
install_prefix:  this directory must end with /cppad_mixed or /dismod_at
build_type:      this must be either debug or release
EOF
	exit 1
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
	-e "s|/cppad_mixed$dollar|/cppad_mixed.$build_type|" \
	-e "s|/dismod_at$dollar|/dismod_at.$build_type|" `
echo "$total_prefix"
if [ "$total_prefix" == "$install_prefix" ]
then
	echo "$program_name: install prefix does not end in /cppad_mixed"
	exit 1
fi
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
