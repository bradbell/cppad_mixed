#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-19 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/check_all.sh" ]
then
	echo 'bin/check_all.sh must be executed from its parent directory'
	exit 1
fi
# ---------------------------------------------------------------------------
if ! grep "build_type='debug'" bin/run_cmake.sh > /dev/null
then
	echo 'bin/check_all.sh: bin/run_cmake.sh build_type not debug'
	exit 1
fi
if ! grep "optimize_cppad_function='no'" bin/run_cmake.sh > /dev/null
then
	echo 'bin/check_all.sh: bin/run_cmake.sh optimize_cppad_function not no'
	exit 1
fi
# ---------------------------------------------------------------------------
release='no'
ldlt_eigen='no'
while [ "$1" != '' ]
do
	if [ "$1" == '--help' ]
	then
		cat << EOF
usage: bin/check_all.sh \\
	[--help] \\
	[--release] \\
	[--ldlt_eigen]
EOF
		exit 0
	fi
	if [ "$1" == '--release' ]
	then
		release='yes'
	elif [ "$1" == '--ldlt_eigen' ]
	then
		ldlt_eigen='yes'
	else
		echo "'$1' is an invalid option"
		bin/check_all.sh --help
		exit 1
	fi
	shift
done
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
# run bin/check_*.sh and ~bradbell/bin/check_copyright.sh
list=`ls bin/check_*.sh`
for script in check_copyright.sh $list
do
	if [ "$script" != 'bin/check_all.sh' ] \
	&& [ "$script" != 'bin/check_install.sh' ]
	then
		$script
	fi
done
#
# check version number
version.sh check
#
# check latex in omhelp
echo_eval run_omhelp.sh -xml dev
echo_eval run_omhelp.sh doc
# -----------------------------------------------------------------------------
# cppad_prefix
cmd=`grep '^cppad_prefix=' bin/run_cmake.sh`
eval $cmd
#
installed_include_dir="$cppad_prefix/include/cppad/mixed"
if [ -e "$installed_include_dir" ]
then
	echo_eval rm -rf $installed_include_dir
fi
# -----------------------------------------------------------------------------
# run_cmake.sh
flags=''
if [ "$release" == 'yes' ]
then
	flags="$flags --release"
fi
if [ "$ldlt_eigen" == 'yes' ]
then
	flags="$flags --ldlt_eigen"
fi
echo "bin/run_cmake.sh $flags >& cmake.log"
bin/run_cmake.sh $flags >& cmake.log
# ----------------------------------------------------------------------------
cd build
#
for target in check speed install
do
	echo "make $target >& $target.log"
	make $target >& ../$target.log
done
cd ..
# -----------------------------------------------------------------------------
check_install_failed='no'
if [ "$release" == 'yes' ]
then
	sed -i bin/run_cmake.sh -e "s|^build_type=.*|build_type='release'|"
	if ! bin/check_install.sh
	then
		check_install_failed='yes'
	fi
	sed -i bin/run_cmake.sh -e "s|^build_type=.*|build_type='debug'|"
else
	if ! bin/check_install.sh
	then
		check_install_failed='yes'
	fi
fi
if [ "$check_install_failed" == 'yes' ]
then
	echo 'bin/check_all.sh: Error'
	exit 1
fi
# -----------------------------------------------------------------------------
for target in cmake check speed install
do
	if grep -i 'warning:' $target.log
	then
		echo "bin/run_check_all.sh: $target.log is has warnings."
		exit 1
	fi
done
# -----------------------------------------------------------------------------
echo 'check_all.sh: OK'
exit 0
