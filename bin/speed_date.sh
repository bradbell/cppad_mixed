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
if [ "$0" != 'bin/speed_date.sh' ]
then
	echo 'bin/speed_date.sh must be run from its parent directory'
	exit 1
fi
if [ "$3" == '' ]
then
	echo 'usage: bin/speed_date.sh branch yyyymmdd quasi_fixed'
	exit 1
fi
branch="$1"
date="$2"
quasi_fixed="$3"
if [ "$quasi_fixed" != 'yes' ] && [ "$quasi_fixed" != 'no' ]
then
	echo 'speed_date.sh: quasi_fixed is not yes or no'
	exit 1
fi
current=`git branch | sed -e '/^[^*]/d' -e 's/^\* *//'`
if [ "$current" != 'master' ]
then
	echo 'bin/speed_branch.sh: must be executed from master branch'
	exit 0
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
	if [ -e "build/$date.$program" ]
	then
		echo_eval rm build/$date.$program
	fi
done
# -----------------------------------------------------------------------------
edit_failed='no'
#
# start with a clean copy of source
git checkout --quiet $branch
hash=`git_date2hash.sh $date`
git checkout $hash
#
# setup to build optimized versions
sed -i bin/run_cmake.sh \
	-e "s|^build_type=.*|build_type='release'|" \
	-e "s|^optimize_cppad_function=.*|optimize_cppad_function='yes'|"
if ! grep "^build_type='release'" bin/run_cmake.sh > /dev/null
then
	edit_failed='build_type'
fi
if ! grep "^optimize_cppad_function='yes'" bin/run_cmake.sh > /dev/null
then
	edit_failed='optimize_cppad_function'
fi
#
bin/run_cmake.sh
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
	if [ "$program" == 'ar1_xam' ]
	then
		sed -i bin/$program.sh \
			-e 's|^number_random=.*|number_random=90000|'
		if ! grep "^number_random=90000" bin/$program.sh > /dev/null
		then
			edit_failed='random_number'
		fi
		#
		sed -i speed/$program.cpp \
			-e 's|\(std::fabs(estimate_ratio\[j\])\) *<.*|\1 < 10.0;|'
		if ! grep "std::fabs(estimate_ratio\[j\]) < 10.0" \
			speed/$program.cpp > /dev/null
		then
			edit_failed='estimate_ratio'
		fi
	fi
	if [ "$edit_failed" != 'no' ]
	then
		echo "bin/speed_date.sh: edit failed: $edit_failed"
		exit 1
	fi
	# -------------------------------------------------------------------
	#
	cd build; make $program; cd ..
	cp build/speed/$program build/$date.$program
	#
	echo "bin/$program.sh normal > build/$date.$program.out"
	bin/$program.sh normal > build/$date.$program.out
done
git reset --hard
# -----------------------------------------------------------------------------
git checkout master
# -----------------------------------------------------------------------------
echo 'bin/speed_date.sh: results are in'
echo "build/ar1_xam.$date1.out,     build/ar1_xam.$date2.out"
echo "build/capture_xam.$date1.out, build/capture_xam.$date2.out"
echo 'speed_date.sh: OK'
exit 0
