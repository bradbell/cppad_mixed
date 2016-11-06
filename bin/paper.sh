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
if [ "$0" != "bin/paper.sh" ]
then
	echo "bin/paper.sh: must be executed from its parent directory"
	exit 1
fi
# ---------------------------------------------------------------------------
check=`bin/version.sh get`
if [ "$check" != "$version" ]
then
	echo 'version and hash_code do not agree'
	exit 1
fi
# ---------------------------------------------------------------------------
# ar1_xam commanmd line options
random_seed='0'
number_random='1000'
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
		file="ar1_xam_${a}_${c}.out"
		use=''
		if [ ! -e "build.release/speed/$file" ]
		then
			use='n'
		else
			while [ "$use" != 'y' ] && [ "$use" != 'n' ]
			do
				read -p "use existing $file [y/n] ?" use
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
			cd build.release
			make speed_ar1_xam
			mv speed/ar1_xam speed/ar1_xam_${a}_${c}
			cd speed
			echo "ar1_xam_${a}_${c} > build.release/speed/$file"
			./ar1_xam_${a}_${c} > $file \
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
		fi
	done
done
# ---------------------------------------------------------------------------
# restore bin/run_camke.sh and bin/check_install.sh
git checkout bin/run_cmake.sh
git checkout bin/check_install.sh
# ---------------------------------------------------------------------------
echo '      (atmoic, checkpoint)'
echo 'name: (yes,yes) (yes,no) (no,yes) (no,no)'
list='
	initialize_bytes
	initialize_seconds
	optimize_fixed_seconds
	information_mat_seconds
'
for name in $list
do
	line="$name"
	for a in yes no
	do
		for c in yes no
		do
			file="build.release/speed/ar1_xam_${a}_${c}.out"
			value=`sed < $file -n -e "/^$name *=/p" | sed  -e 's|.*= *||'`
			line="$line $value"
		done
	done
	echo "$line"
done
# ---------------------------------------------------------------------------
echo 'compare.sh: OK'
