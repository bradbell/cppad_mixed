#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
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
branch=`git branch | grep '^\*' | sed -e 's|^\* *||'`
if [ "$branch" != 'master' ]
then
	echo 'bin/speed_branch.sh must be run from master branch'
	exit 1
fi
if [ "$3" == '' ]
then
	echo 'usage: bin/speed_branch.sh branch1 branch2 quasi_fixed'
	exit 1
fi
branch1="$1"
branch2="$2"
quasi_fixed="$3"
if [ "$quasi_fixed" != 'yes' ] && [ "$quasi_fixed" != 'no' ]
then
	echo 'speed_branch.sh: quasi_fixed is not yes or no'
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
	diff=`git diff $branch1 $branch2 -- bin/$program`
	if [ "$diff" != '' ]
	then
		echo "bin/speed_branch.sh: bin/$program.sh has changed"
		exit 1
	fi
	for ext in $branch1 $branch2
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
for branch in $branch1 $branch2
do
	git checkout $branch
	bin/run_cmake.sh --release
	for program in ar1_xam capture_xam
	do
		cd build; make $program; cd ..
		cp build/speed/$program build/$branch.$program
		#
		sed -i bin/$program.sh \
			-e 's|^random_seed=.*|random_seed=123|' \
			-e "s|^quasi_fixed=.*|quasi_fixed=$quasi_fixed|"
		#
		if [ "$program" == 'ar1_xam' ]
		then
			sed -i bin/$program.sh \
				-e 's|^number_random=.*|number_random=90000|'
		fi
		#
		echo "bin/$program.sh normal > build/$branch.$program.out"
		bin/$program.sh normal > build/$branch.$program.out
	done
	git reset --hard
done
git checkout master
# -----------------------------------------------------------------------------
echo 'bin/speed_branch.sh: results are in'
echo "build/ar1_xam.$branch1.out,     build/ar1_xam.$branch2.out"
echo "build/capture_xam.$branch1.out, build/capture_xam.$branch2.out"
echo 'speed_branch.sh: OK'
exit 0
