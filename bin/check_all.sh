#! /bin/bash -e
# $Id:$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != "bin/check_all.sh" ]
then
	echo 'bin/check_all.sh [d|r]'
	echo 'must be executed from its parent directory'
	exit 1
fi
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
input="$1"
while [ "$input" != 'd' ]  &&
      [ "$input" != 'r' ]  &&
      [ "$input" != 'da' ] &&
      [ "$input" != 'ra' ]
do
	msg='debug [d], release [r], debug & atomic_chol [da]'
	msg="$msg, release & atomic_chol [ra] ?"
	read -p "$msg" input
done
# -----------------------------------------------------------------------------
list=`ls bin/check_*.sh`
for script in $list
do
	if [ "$script" != 'bin/check_all.sh' ] \
	&& [ "$script" != 'bin/check_install.sh' ]
	then
		$script
	fi
done
# ----------------------------------------------------------------------------
if [ "$input" == 'd' ]
then
	build_type='debug'
	bin/run_cmake.sh
elif [ "$input" == 'da' ]
then
	build_type='debug'
	bin/run_cmake.sh --use_sparse_ad_cholesky
elif [ "$input" == 'r' ]
then
	build_type='release'
	bin/run_cmake.sh --release
elif [ "$input" == 'ra' ]
then
	build_type='release'
	bin/run_cmake.sh --release --use_sparse_ad_cholesky
else
	echo 'error in check_all.sh script'
	exit 1
fi
# ----------------------------------------------------------------------------
bin/run_omhelp.sh xml
#
sed -i bin/check_install.sh \
	-e "s|^build_type=.*|build_type='$build_type'|"
#
cd build
make check
make speed
make install
cd ..
bin/check_install.sh
#
if [ "$build_type" == 'release' ]
then
	sed -i bin/check_install.sh \
		-e "s|^build_type=.*|build_type='debug'|"
fi
#
echo 'check_all.sh: OK'
