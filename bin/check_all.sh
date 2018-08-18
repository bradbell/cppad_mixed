#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
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
build_type="$1"
if [ "$build_type" != 'debug' ]  && [ "$build_type" != 'release' ]
then
	echo 'bin/check_all.sh (debug|release)'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
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
echo_eval run_omhelp.sh -xml doc
# -----------------------------------------------------------------------------
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
if [ "$build_type" == 'release' ]
then
	echo 'bin/run_cmake.sh --release --optimize_cppad_function >& cmake.log'
	bin/run_cmake.sh --release --optimize_cppad_function >& cmake.log
else
	echo 'bin/run_cmake.sh >& cmake.log'
	bin/run_cmake.sh >& cmake.log
fi
#
cd build
#
for target in check speed install
do
	echo "make $target >& $target.log"
	make $target >& ../$target.log
done
cd ..
bin/check_install.sh
# -----------------------------------------------------------------------------
for target in cmake check speed install
do
	if grep -i 'warningL:' $target.log
	then
		echo "bin/run_check_all.sh: $target.log is has warnings."
		exit 1
	fi
done
# -----------------------------------------------------------------------------
echo 'check_all.sh: OK'
exit 0
