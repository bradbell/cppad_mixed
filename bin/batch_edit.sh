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
'
#
rename_cmd='s|sparse_ad_cholesky|sparse_ad_chol|'
#
cat << EOF > junk.sed
/^\\t*mixed_derived(\$/! b mixed_object_on_multiple_lines
# define derived class constructor in defined here
: loop
N
/cppad_mixed(/! b loop
s|\\n\\t*const *sparse_mat_info& *A_info *,||
s| *,\\n\\t*const *sparse_mat_info& *A_info||
s|,\\n\\t*CppAD::mixed::sparse_mat_info& *A_info *||
s|cppad_mixed(\\([^,]*,[^,]*,[^,]*\\), *A_info *)|cppad_mixed(\\1)|
s|cppad_mixed(1, 0, false, A_info)|cppad_mixed(1, 0, false)|
b done
# ---------------------------------------------------------------------------
: mixed_object_on_multiple_lines
/^\\t*mixed_derived *mixed_object( *\$/! b normal
N
N
s|mixed_object(\\([^,]*,[^,]*\\), *A_info *\\([),]\\)|mixed_object(\\1\\2|
s|mixed_object(\\([^,]*,[^,]*,[^,]*\\), *A_info *\([),]\\)|mixed_object(\\1\\2|
s|mixed_object(\\([^,]*,[^,]*,[^,]*,[^,]*\\), *A_info *\([),]\\)|mixed_object(\\1\\2|
b done
# ----------------------------------------------------------------------------
: normal
s|/derived_ctor/A_info/|/initialize/A_info/|
s|n_random, quasi_fixed, A_info|n_random, quasi_fixed|
#
s|initialize(fixed_in, random_in)|initialize(fixed_in, random_in, A_info)|
s|initialize(theta_in, u_in)|initialize(theta_in, u_in, A_info)|
s|initialize(fixed_vec, random_vec)|initialize(fixed_vec, random_vec, A_info)|
s|initialize(fixed_vec, random_in)|initialize(fixed_vec, random_in, A_info)|
s|initialize(theta, u)|initialize(theta, u, A_info)|
#
s|mixed_object(\\([^,]*,[^,]*\\), *A_info *\\([),]\\)|mixed_object(\\1\\2|
s|mixed_object(\\([^,]*,[^,]*,[^,]*\\), *A_info *\([),]\\)|mixed_object(\\1\\2|
s|mixed_object(\\([^,]*,[^,]*,[^,]*,[^,]*\\), *A_info *\([),]\\)|mixed_object(\\1\\2|
s|^\\(\\t*n_fixed *, *n_random *,\\) *A_info *,|\\1|
: done
# ----------------------------------------------------------------------------
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
