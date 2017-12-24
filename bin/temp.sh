#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-17 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
if [ "$0" != 'bin/temp.sh' ]
then
	echo 'bin/temp.sh must be run from its parent directory'
	exit 1
fi
if [ "$1" == '' ]
then
	echo 'usage: bin/temp.sh file'
	exit 1
fi
file="$1"
# -----------------------------------------------------------------------------
# bash function that echos and executes a command
echo_eval() {
	echo $*
	eval $*
}
# -----------------------------------------------------------------------------
cat << EOF > junk.sed
s|Copyright (C) 2014-1\\([4-6]\\)|Copyright (C) 2014-17|
# ----------------------------------------------------------------------------
/^\\t\\tvirtual a1_vector fix_likelihood(/! b two
: one
N
/\\n\\t\\t}/! b one
s|\\t\\tvirtual a1_vector fix_likelihood(|\\t\\ttemplate <typename Vector>\\
\\t\\tVector template_fix_likelihood(|
#
s|\\n\\t\\t{|&\\ttypedef typename Vector::value_type scalar;\\n\\n\\t\\t|
s|a1_vector|Vector|g
s|a1_double|scalar|g
#
s|\\n\\t\\t}|&\\
\\t\\t// a1_vector version of fix_likelihood\\
\\t\\tvirtual a1_vector fix_likelihood(const a1_vector\\& fixed_vec)\\
\\t\\t{	return template_fix_likelihood( fixed_vec ); }\\
\\t\\t// Delete this a2_vector version before commiting\\
\\t\\tvirtual a2_vector fix_likelihood(const a2_vector\\& fixed_vec)\\
\\t\\t{	return template_fix_likelihood( fixed_vec ); }|
# ----------------------------------------------------------------------------
: two
/^\\t\\tvirtual a2_vector ran_likelihood(/! b four
:three
N
/\\n\\t\\t}/! b three
s|\\t\\tvirtual a2_vector ran_likelihood(|\\t\\ttemplate <typename Vector>\\
\\t\\tVector template_ran_likelihood(|
#
s|\\n\\t\\t{|&\\ttypedef typename Vector::value_type scalar;\\n\\n\\t\\t|
s|a2_vector|Vector|g
s|a2_double|scalar|g
#
s|\\n\\t\\t}|&\\
\\t\\t// a2_vector version of ran_likelihood\\
\\t\\tvirtual a2_vector ran_likelihood(\\
\\t\\t	const a2_vector\\& fixed_vec, const a2_vector\\& random_vec\\
\\t\\t)\\
\\t\\t{	return template_ran_likelihood( fixed_vec, random_vec ); }\\
\\t\\t// Delete this a3_vector before commiting\\
\\t\\tvirtual a3_vector ran_likelihood(\\
\\t\\t	const a3_vector\\& fixed_vec, const a3_vector\\& random_vec\\
\\t\\t)\\
\\t\\t{	return template_ran_likelihood( fixed_vec, random_vec ); }|
# ----------------------------------------------------------------------------
: four
/^\\t\\tvirtual a1_vector fix_constraint(/! b six
: five
N
/\\n\\t\\t}/! b five
s|\\t\\tvirtual a1_vector fix_constraint(|\\t\\ttemplate <typename Vector>\\
\\t\\tVector template_fix_constraint(|
#
s|\\n\\t\\t{|&\\ttypedef typename Vector::value_type scalar;\\n\\n\\t\\t|
s|a1_vector|Vector|g
s|a1_double|scalar|g
#
s|\\n\\t\\t}|&\\
\\t\\t// a1_vector version of fix_constraint\\
\\t\\tvirtual a1_vector fix_constraint(const a1_vector\\& fixed_vec)\\
\\t\\t{	return template_fix_constraint( fixed_vec ); }\\
\\t\\t// Delete this a2_vector version before commiting\\
\\t\\tvirtual a2_vector fix_constraint(const a2_vector\\& fixed_vec)\\
\\t\\t{	return template_fix_constraint( fixed_vec ); }|
# ----------------------------------------------------------------------------
: six
EOF
# -----------------------------------------------------------------------------
git checkout $file
echo_eval sed -i $file -f junk.sed