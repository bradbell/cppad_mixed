#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-20 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/speed_test.sh' ]
then
	echo 'bin/speed_test.sh must be run from its parent directory'
	exit 1
fi
branch=`git branch | grep '^\*' | sed -e 's|^\* *||'`
if [ "$branch" != 'master' ]
then
	echo 'bin/speed_test.sh must be run from master branch'
	exit 1
fi
if [ "$3" == '' ]
then
	echo 'usage: bin/speed_test.sh branch1 (branch2|new) quasi_fixed'
	echo 'where branch1 and branch2 are git branches and new means'
	echo 'branch1 with the changes corresponding to get_new.sh from'
	exit 1
fi
branch1="$1"
branch2="$2"
quasi_fixed="$3"
if [ "$quasi_fixed" != 'yes' ] && [ "$quasi_fixed" != 'no' ]
then
	echo 'speed_test.sh: quasi_fixed is not yes or no'
	exit 1
fi
if [ "$branch1" == 'new' ]
then
	echo 'speed_test.sh: branch1 cannot be "new"'
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
# erase old version of programs
for program in ar1_xam capture_xam
do
	for ext in $branch1 $branch2
	do
		if [ -e "build.release/$program.$ext" ]
		then
			echo_eval rm build.release/$program.$ext
		fi
		if [ -e "build.release/$ext.$program" ]
		then
			echo_eval rm build.release/$ext.$program
		fi
	done
done
# -----------------------------------------------------------------------------
for branch in $branch1 $branch2
do
	edit_failed='no'
	#
	# start with a clean copy of branch source
	if [ $branch == 'new' ]
	then
		git_new.sh from
	else
		git checkout $branch
	fi
	#
	echo_eval bin/run_cmake.sh --release --optimize_cppad_function
	#
	# for each test
	for program in ar1_xam capture_xam
	do
		# -------------------------------------------------------------------
		# changes to source code
		sed -i bin/$program.sh \
			-e 's|^random_seed=.*|random_seed=123|' \
			-e "s|^quasi_fixed=.*|quasi_fixed=$quasi_fixed|"
		if ! grep "^random_seed=123" bin/$program.sh > /dev/null
		then
			edit_failed='random_seed'
		fi
		if ! grep "^quasi_fixed=$quasi_fixed" bin/$program.sh > /dev/null
		then
			edit_failed='quasi_fixed'
		fi
		#
		if [ "$edit_failed" != 'no' ]
		then
			echo "bin/speed_test.sh: edit failed: $edit_failed"
			exit 1
		fi
		# -------------------------------------------------------------------
		#
		cd build; make $program; cd ..
		cp build/speed/$program build/$branch.$program
		#
		echo "bin/$program.sh normal > build/$branch.$program.out"
		bin/$program.sh normal > build/$branch.$program.out
	done
	git reset --hard
done
git checkout master
# -----------------------------------------------------------------------------
echo 'bin/speed_test.sh: results are in:'
for branch in $branch1 $branch2
do
	for program in ar1_xam capture_xam
	do
		echo "	build/$branch.$program.out"
	done
done
echo 'speed_test.sh: OK'
exit 0
