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
ok='yes'
if [ "$1" != 'normal' ] && [ "$1" != 'massif' ] && [ "$1" != 'callgrind' ]
then
	ok='no'
fi
if [ "$2" != 'yes' ] && [ "$2" != 'no' ]
then
	ok='no'
fi
if [ "$ok" == 'no' ]
then
	echo 'usage: bin/speed_new.sh test2run quasi_fixed'
	echo 'test2run:    is normal, massif, or callgrid'
	echo 'quasi_fixed: is yes or no'
	exit 1
fi
test2run="$1"
quasi_fixed="$2"
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
bin/run_cmake.sh --release --optimize_cppad_function
# -----------------------------------------------------------------------------
for program in ar1_xam capture_xam
do
	diff=`git diff bin/$program.sh`
	if [ "$diff" != '' ]
	then
		echo "bin/speed_new.sh: bin/$program.sh has changed"
		exit 1
	fi
	for ver in old new
	do
		if [ -e "build/$ver.$program" ]
		then
			echo_eval rm build/$ver.$program
		fi
	done
done
# -----------------------------------------------------------------------------
# old version of code
for ver in old new
do
	for program in ar1_xam capture_xam
	do
		# make sure quasi_fixed is as expected
		sed -i bin/$program.sh \
			-e "s|^quasi_fixed=.*|quasi_fixed=$quasi_fixed|"
		#
		cd build; make $program; cd ..
		cp build/speed/$program build/$ver.$program
		#
		echo "bin/$program.sh $test2run > build/$ver.$program.out"
		bin/$program.sh $test2run > build/$ver.$program.out
		if [ "$test2run" == 'massif' ]
		then
			mv build/speed/massif.out build/$ver.massif.out
			sed -i build/$ver.$program.out \
				-e 's|build/speed/massif.out|build/massif.out|'
		fi
		git checkout bin/$program.sh
	done
	# This may overwrite bin/$program.sh
	git_new.sh from
done
# -----------------------------------------------------------------------------
echo 'speed_new.sh: results are in'
echo 'build/old.ar1_xam.out,     build/new.ar1_xam.out'
echo 'build/old.capture_xam.out, build/new.capture_xam.out'
echo 'speed_new.sh: OK'
exit 0
