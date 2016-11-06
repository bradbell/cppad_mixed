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
version='20161104'
hash_code='47f22afd7d2495b888404922ec3fd0eedb0f277e'
if [ "$0" != "bin/paper.sh" ]
then
	echo "bin/paper.sh: must be executed from its parent directory"
	exit 1
fi
# ---------------------------------------------------------------------------
git checkout --quiet $hash_code
check=`bin/version.sh get`
if [ "$check" != "$version" ]
then
	echo 'version and hash_code do not agree'
	exit 1
fi
# ---------------------------------------------------------------------------
# ar1_xam commanmd line options
random_seed='0'
number_random='500'
quasi_fixed='false'
trace_optimize_fixed='false'
ipopt_solve='false'
bool_sparsity='false'
hold_memory='false'
derivative_test='false'
start_near_solution='false'
# ---------------------------------------------------------------------------
# modify bin/run_cmake.sh and build different versions of ar1_xam
git checkout bin/run_cmake.sh
sed -i bin/run_cmake.sh \
	-e "s|^\(build_type\)=.*|\1='release'|" \
	-e "s|^\(optimize_cppad_function\)=.*|\1='YES'|"
# ---------------------------------------------------------------------------
for atomic in YES NO
do
	for checkpoint in YES NO
	do
		a=`echo "$atomic" | tr A-Z a-z`
		c=`echo "$checkpoint" | tr A-Z a-z`
		use=''
		if [ ! -e "build/speed/ar1_xam_${a}_${c}" ]
		then
			use='n'
		else
			while [ "$use" != 'y' ] && [ "$use" != 'n' ]
			do
				read -p "use existing ar1_xam_${a}_${c} [y/n] ?" use
			done
		fi
		if [ "$use" == 'n' ]
		then
			sed -i bin/run_cmake.sh \
				-e "s|^\(use_atomic_cholesky\)=.*|\1='$atomic'|" \
				-e "s|^\(checkpoint_newton_step\)=.*|\1='$checkpoint'|"
			bin/run_cmake.sh > paper.tmp
			#
			if ! grep "use_atomic_cholesky = $atomic" paper.tmp
			then
				echo 'error in compare.sh'
				exit 1
			fi
			#
			if ! grep "checkpoint_newton_step = $checkpoint" paper.tmp
			then
				echo 'error in compare.sh'
				exit 1
			fi
			rm paper.tmp
			#
			cd build
			make speed_ar1_xam
			mv speed/ar1_xam speed/ar1_xam_${a}_${c}
			cd ..
		fi
		#
		cd build/speed
		echo "ar1_xam_${a}_${c} > build/speed/ar1_xam_${a}_${c}.out"
		./ar1_xam_${a}_${c} > ar1_xam_${a}_${c}.out \
			$random_seed \
			$number_random \
			$quasi_fixed \
			$trace_optimize_fixed \
			$ipopt_solve \
			$bool_sparsity \
			$hold_memory \
			$derivative_test \
			$start_near_solution
		cd ../..
	done
done
# ---------------------------------------------------------------------------
# restore bin/run_camke.sh and bin/check_install.sh
git checkout bin/run_cmake.sh
git checkout bin/check_install.sh
# ---------------------------------------------------------------------------
echo 'compare.sh: OK'
