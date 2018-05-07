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
# This scipt demonstrates that changing cppad-20171215 to cppad-20180503
# uses less mmeory and does not affect speed.
declare -A branch
branch[one]='6894f1'
branch[two]='master'
# -----------------------------------------------------------------------------
if [ "$0" != 'bin/speed_temp.sh' ]
then
	echo 'bin/speed_temp.sh must be run from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
build_dir='build/speed'
for version in one two
do
	rebuild='no'
	for program in capture_xam ar1_xam
	do
		if [ ! -e $build_dir/$program.$version ]
		then
			rebuild='yes'
		fi
	done
	if [ "$rebuild" == 'yes' ]
	then
		echo_eval git checkout ${branch[$version]}
		#
		# temporary change to bin/run_cmake.sh
		sed -i bin/run_cmake.sh \
			-e "s|^build_type=.*|build_type='release'|"
		#
		echo_eval bin/install_cppad.sh
		echo_eval bin/run_cmake.sh --release
		#
		echo_eval pushd $build_dir
		for program in capture_xam ar1_xam
		do
			echo_eval make $program
			echo_eval mv $program $program.$version
		done
		popd
		#
		# undo changes in bin/run_cmake.sh
		git checkout bin/run_cmake.sh
	fi
	#
	for program in capture_xam ar1_xam
	do
		# temporary change to bin/$program.sh
		git show master:bin/$program.sh | sed \
			-e "s|random_seed='0'|random_seed='123'|" \
			-e "s|number_locations='30'|number_locations='30'|" \
			-e "s|number_random='10000'|number_random='10000'|" \
			-e "s|\.\\/$program |./$program.$version |" \
			> bin/$program.sh
		#
		# create a phony $build_dir/$program to avoid error in bin/$program.sh
		# and to make sure it does not get executed
		if [ -e "$build_dir/$program" ]
		then
			rm $build_dir/$program
		fi
		touch $build_dir/$program
		#
		echo "bin/$program.sh normal > $build_dir/$program.$version.out"
		if bin/$program.sh normal > $build_dir/$program.$version.out
		then
			echo "$program: Test passed"
		else
			echo "$program: Test failed"
		fi
		#
		# undo changes to bin/$program.sh
		git checkout bin/$program.sh
	done
	#
	pushd $build_dir
	for program in capture_xam ar1_xam
	do
		cmd=`grep "^\\.\\/$program.$version " $program.$version.out`
		echo "massif.sh $cmd >& /dev/null"
		if massif.sh $cmd >& /dev/null
		then
			echo "$program: Test passed"
		else
			echo "$program: Test failed"
		fi
		echo_eval mv massif.out $program.$version.massif
	done
	popd
done
echo "one = ${branch[one]} two = ${branch[two]}"
for program in capture_xam ar1_xam
do
	for version in one two
	do
		seconds=`grep optimize_fixed_seconds $build_dir/$program.$version.out`
		echo "$program $version: $seconds"
	done
done
