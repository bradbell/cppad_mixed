#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/speed_branch.sh' ]
then
	echo 'bin/speed_branch.sh must be run from its parent directory'
	exit 1
fi
if [ "$1" == '' ] || [ "$2" != '' ]
then
	echo 'usage: bin/speed_branch.sh branch_or_hash'
	exit 1
fi
branch="$1"
if [ "$branch" == 'master' ]
then
	echo 'speed_branch.sh: branch_or_hash cannot be master'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
# set random seed to 123 so same for old and new
for program in ar1_xam capture_xam
do
	diff=`git diff master $branch -- bin/$program`
	if [ "$diff" != '' ]
	then
		echo "bin/speed_branch.sh: bin/$program.sh has changed"
		exit 1
	fi
	for ext in master $branch
	do
		if [ -e "build/$program.$ext" ]
		then
			echo_eval rm build/$program.$ext
		fi
		if [ -e "build/$ext.$program" ]
		then
			echo_eval rm build/$ext.$program
		fi
	done
done
# -----------------------------------------------------------------------------
# branch version of code
git checkout $branch
bin/run_cmake.sh --optimize_cppad_function --release
for program in ar1_xam capture_xam
do
	cd build; make $program; cd ..
	cp build/speed/$program build/$branch.$program
	#
	sed -i bin/$program.sh -e 's|^random_seed=.*|random_seed=123|'
	echo "bin/$program.sh normal > build/$program.$branch"
	bin/$program.sh normal > build/$program.$branch
done
git reset --hard
# -----------------------------------------------------------------------------
# master version of code
git checkout master
bin/run_cmake.sh --optimize_cppad_function --release
for program in ar1_xam capture_xam
do
	cd build; make $program; cd ..
	cp build/speed/$program build/master.$program
	#
	sed -i bin/$program.sh -e 's|^random_seed=.*|random_seed=123|'
	echo "bin/$program.sh normal > build/$program.master"
	bin/$program.sh normal > build/$program.master
done
git reset --hard
# -----------------------------------------------------------------------------
echo 'bin/speed_branch.sh: results are in'
echo "build/ar1_xam.master,     build/ar1_xam.$branch"
echo "build/capture_xam.master, build/capture_xam.$branch"
echo 'speed_branch.sh: OK'
exit 0
