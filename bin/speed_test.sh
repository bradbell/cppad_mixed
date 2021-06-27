#! /bin/bash -e
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-21 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#        GNU Affero General Public License version 3.0 or later
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
if [ "$2" == '' ]
then
	echo 'usage: bin/speed_test.sh (branch|new) quasi_fixed'
	exit 1
fi
branch="$1"
quasi_fixed="$2"
if [ "$quasi_fixed" != 'yes' ] && [ "$quasi_fixed" != 'no' ]
then
	echo 'speed_test.sh: quasi_fixed is not yes or no'
	exit 1
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
	if [ -e "build.release/$program.$branch" ]
	then
		echo_eval rm build.release/$program.$branch
	fi
	if [ -e "build.release/$branch.$program" ]
	then
		echo_eval rm build.release/$branch.$program
	fi
done
# -----------------------------------------------------------------------------
edit_failed='no'
#
# start with a clean copy of branch source
git reset --hard
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
git checkout master
# -----------------------------------------------------------------------------
echo 'bin/speed_test.sh: results are in:'
for program in ar1_xam capture_xam
do
	echo "	build/$branch.$program.out"
done
echo 'speed_test.sh: OK'
exit 0
