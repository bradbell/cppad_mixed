#! /bin/bash -e
# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
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
for program in ar1_xam capture_xam
do
	sed -i bin/$program.sh -e 's|^random_seed=.*|random_seed=123|'
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
	echo "bin/$program.sh normal > build/$program.old"
	bin/$program.sh normal > build/$program.old
done
# -----------------------------------------------------------------------------
# new version of code
git_new.sh from
for program in ar1_xam capture_xam
do
	cd build; make $program; cd ..
	echo "bin/$program.sh normal > build/$program.new"
	bin/$program.sh normal > build/$program.new
done
# -----------------------------------------------------------------------------
echo 'speed_new.sh: results are in'
echo 'build/ar1_xam.old,     build/ar1_xam.new'
echo 'build/capture_xam.old, build/capture_xam.new'
echo 'speed_new.sh: OK'
exit 0

