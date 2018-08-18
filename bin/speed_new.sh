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
if [ "$0" != 'bin/speed_new.sh' ]
then
	echo 'bin/speed_new.sh must be run from its parent directory'
	exit 1
fi
test2run='normal'
quasi_fixed='yes'
while [ "$1" != '' ]
do
	if [ "$1" == '--callgrind' ]
	then
		if [ "$test2run" == 'massif' ]
		then
			echo 'usage: speed_new.sh: cant specify both callgrind and massif'
			exit 1
		fi
		test2run='callgrind'
	elif [ "$1" == '--massif' ]
	then
		if [ "$test2run" == 'callgrind' ]
		then
			echo 'usage: speed_new.sh: cant specify both callgrind and massif'
			exit 1
		fi
		test2run='massif'
	elif [ "$1" == '--newton' ]
	then
		quasi_fixed='no'
	else
		echo 'usage: bin/speed_new.sh: [--callgrind] [--massif] [--newton]'
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
if [ ! -e new ]
then
	echo 'speed_new.sh: The directory ./new does not exist'
	exit 1
fi
# -----------------------------------------------------------------------------
git reset --hard
bin/run_cmake.sh --release
# -----------------------------------------------------------------------------
# set random seed to 123 so same for old and new
for program in ar1_xam capture_xam
do
	diff=`git diff bin/$program.sh`
	if [ "$diff" != '' ]
	then
		echo "bin/speed_new.sh: bin/$program.sh has changed"
		exit 1
	fi
	git checkout bin/$program.sh
	sed -i bin/$program.sh -e "s|^quasi_fixed=.*|quasi_fixed=$quasi_fixed|"
	for ext in old new
	do
		if [ -e "build/$program.$ext" ]
		then
			echo_eval rm build/$program.$ext
		fi
	done
done
# -----------------------------------------------------------------------------
# old version of code
for program in ar1_xam capture_xam
do
	cd build; make $program; cd ..
	echo "bin/$program.sh $test2run > build/$program.old"
	bin/$program.sh $test2run > build/$program.old
	if [ "$test2run" == 'massif' ]
	then
		mv build/speed/massif.out build/massif.old
		sed -i build/$program.old \
			-e 's|build/speed/massif.out|build/massif.old|'
	fi
done
# -----------------------------------------------------------------------------
# new version of code
git_new.sh from
for program in ar1_xam capture_xam
do
	cd build; make $program; cd ..
	echo "bin/$program.sh $test2run > build/$program.new"
	bin/$program.sh $test2run > build/$program.new
	if [ "$test2run" == 'massif' ]
	then
		mv build/speed/massif.out build/massif.new
		sed -i build/$program.new \
			-e 's|build/speed/massif.out|build/massif.new|'
	fi
done
# -----------------------------------------------------------------------------
# restore random seed
for program in ar1_xam capture_xam
do
	git checkout bin/$program.sh
done
# -----------------------------------------------------------------------------
echo 'speed_new.sh: results are in'
echo 'build/ar1_xam.old,     build/ar1_xam.new'
echo 'build/capture_xam.old, build/capture_xam.new'
echo 'speed_new.sh: OK'
exit 0
