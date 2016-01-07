$Id:$
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
/*
$begin ran_con_matrix_$$
$spell
$$

$section Random Constraint Matix$$

$head Purpose$$
This variable should be a $cref private$$ member variable.
It is instead a static variable in the $code CppAD::mixed$$ namespace
so the warnings that Eigen generates
do not need to be suppressed by all the routines that include
$code cppad_mixed/cppad_mixed.hpp$$.

$head Declaration$$
This variable has the following declaration
$codep */
	namespace CppAD { namespace mixed {
		Eigen::sparseMatrix<double, Eigen::ColMajor> ran_con_matrix_;
	} }
/*
$end
------------------------------------------------------------------------------
$begin set_ran_con_matrix$$
$spell
$$

$section Setting Random Constraints Matrix$$

$head Syntax$$
$code%CppAD::mixed::set_ran_con_matrix(%n_random%, %row%, %col%, %val%)
%$$

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is the number of
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ in the model.
In the case where there are
$cref/no random effects/cppad_mixed/Problem/No Random Effects/$$,
$icode%n_random% = 0%$$.


$subhead row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
It contains one element each non-zero entry in
$cref/random constraint matrix
	/ran_constraint
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$.
We use the notation $icode%K% = %row%.size()%$$.
For $icode%k% = 0 , %...%, %K%$$,
$icode%row%[%k%]%$$ is the row index of a possibly non-zero entry
in $icode A$$.

$subhead col$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %col%
%$$
and has the same size as $icode row$$.
For $icode%k% = 0 , %...%, %K%$$,
$codei%col%[%k%]%$$ is the column index of a non-zero entry
in $icode A$$.
Note that for each $icode k$$, $codei%col%[%k%]% < %n_rancom%$$.


$subhead val$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val%
%$$
and has the same size as $icode row$$.
For $icode%k% = 0 , %...%, %K%$$,
$icode%val%[%k%]%$$ is the value of a non-zero entry
in $icode A$$.
Note that for each $icode k$$, $codei%val%[%k%]% != 0.0$$.

$subhead ran_con_matrix_$$
The input value of this matrx does not matter.
Upon return it is an $code Eigen$$ sparse matrix representation
of the random constraint matrix $latex A$$.

$end
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void set_ran_con_matrix(
	size_t                       n_random ,
	const CppAD::vector<size_t>& row      ,
	const CppAD::vector<size_t>& col      ,
	const CppAD::vector<size_t>& val      )
{	assert( row.size() == col.size() );
	assert( row.size() == val.size() );

	// number of possibly non-zero entries
	size_t K = row.size();

	// determine maximum row and column index
	size_t max_row_p1 = 0;
	size_t max_col_p1 = 0;
	for(size_t k = 0; k < K; k++)
	{	assert( val[k] != 0.0 );
		max_row_p1 = std::max(max_row_p1, row[k] + 1);
		max_col_p1 = std::max(max_col, col[k] + 1);
	}
	assert( max_col_p1 <= n_random );

	// size of the matrix A
	nrow = max_row_p1;
	ncol = n_random;
	ran_con_matrix_.resize(nrow, ncol);

	// put the values in A
	for(size_t k = 0; k < K; k++)
		ran_con_matrix_.insert(row[k], col[k]) = val[k];

	return;
}
