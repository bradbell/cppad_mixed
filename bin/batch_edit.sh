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
new_directories='
'
rename_files='
'
spell_files='
'
no_change_files='
	test_more/no_random_info.cpp
'
#
rename_cmd='s|private/box_newton_xam.cpp|user/box_newton_xam.cpp|'
#
cat << EOF > junk.sed
s|// Float sqrt_2pi = Float(\\(.*\\) );|// sqrt_2pi =\\1;|
s|^\\tusing CppAD::mixed::sparse_mat_info;|&\\
\\t//\\
\\ttypedef AD<double>    a1_double;\\
\\ttypedef AD<a1_double> a2_double;|
#
/^\\tprivate:/! b one
N
/\/\/ implementation/! b one
s|^\\tprivate:\\n||
# -------------------------------------------------------
: one
#
/^\\t\\ttemplate <class Float>/! b five
N
/vector<Float> implement/! b five
# -------------------------------------------------------
/implement_fix_likelihood/! b three
: two
N
/\\n\\t\\t}\$/! b two
s|^\\t\\ttemplate <class Float>\\n||
s|Float|a1_double|g
s|vector<a1_double> implement_fix_likelihood|virtual &|
s|implement_fix_likelihood|fix_likelihood|
b end
# -------------------------------------------------------
: three
/implement_ran_likelihood/! b end
: four
N
/\\n\\t\\t}\$/! b four
s|^\\t\\ttemplate <class Float>\\n||
s|Float|a2_double|g
s|vector<a2_double> implement_ran_likelihood|virtual &|
s|implement_ran_likelihood|ran_likelihood|
b end
# -------------------------------------------------------
: five
/^\\t\\tvirtual/! b end
: six
N
/return implement/! b six
d
# -------------------------------------------------------
: end
EOF
# -----------------------------------------------------------------------------
if [ "$0" != "bin/batch_edit.sh" ]
then
	echo "bin/batch_edit.sh: must be executed from its parent directory"
	exit 1
fi
# -----------------------------------------------------------------------------
cp bin/batch_edit.sh $HOME/trash/batch_edit.sh
git reset --hard
# ----------------------------------------------------------------------------
# new directories
for dir in $new_directories
do
	if [ ! -e $dir ]
	then
		mkdir $dir
	fi
done
# ----------------------------------------------------------------------------
# rename files
for old in $rename_files
do
	new=`echo $old | sed -e "$rename_cmd"`
	echo_eval git mv $old $new
	echo "s|$old|$new|g" >> junk.sed
done
# ----------------------------------------------------------------------------
list_all=`bin/ls_files.sh`
for file in $list_all
do
	if [ "$file" != 'bin/batch_edit.sh' ]
	then
		sed -f junk.sed $file > junk.$$
		if diff $file junk.$$ > /dev/null
		then
			rm junk.$$
		else
			echo "sed -f junk.sed -i $file"
			if [ -x "$file" ]
			then
				mv junk.$$ $file
				chmod +x $file
			else
				mv junk.$$ $file
			fi
		fi
	fi
done
# ----------------------------------------------------------------------------
# files that have spelling changes
for file in $spell_files
do
	sed -e "$spell_cmd" -i $file
done
# ----------------------------------------------------------------------------
# files that should not change at all
for file in $no_change_files
do
	echo_eval git checkout $file
done
# ----------------------------------------------------------------------------
# files that have hand edited and stored using 'git_new.sh to'
git_new='yes'
for name in add_list mod_list rm_list
do
	if [ ! -e new/$name ]
	then
		git_new='no'
	fi
done
if [ "$git_new" == 'yes' ]
then
	git_new.sh from
fi
# ----------------------------------------------------------------------------
echo_eval cp $HOME/trash/batch_edit.sh bin/batch_edit.sh
echo 'batch_edit.sh: OK'
exit 0
