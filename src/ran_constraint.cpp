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
$begin cppad_mixed_static$$
$spell
$$

$section Declare Some CppAD::mixed Static Variables$$

$head Purpose$$
These variables should be a $cref private$$ member variable.
It is instead a static variable in the $code CppAD::mixed$$ namespace
so the warnings that Eigen generates
do not need to be suppressed by all the routines that include
$code cppad_mixed/cppad_mixed.hpp$$.

$head ran_con_matrix_$$
This is the $cref/random constraint matrix
	/ran_constraint
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$.
It has the following declaration
$codep */
	namespace CppAD { namespace mixed {
		Eigen::sparseMatrix<double, Eigen::ColMajor> ran_con_matrix_;
	} }

/* $$

$head hes_cross_matrix$$
This is the current value of the Hessian cross terms.
$latex f_{u, \theta } ( \theta , u)$$.
$codep */
	namespace CppAD { namespace mixed {
		Eigen::sparseMatrix<double, Eigen::ColMajor> hes_cross_matrix_;
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
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$.
It is assumed that $icode%n_random% > 0%$$.


$subhead row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
It contains one element each possibly non-zero entry in
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
$icode%val%[%k%]%$$ is the value of a possibly non-zero entry
in $icode A$$.

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
{	assert( n_random > 0 );
	assert( row.size() == col.size() );
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
/*
$end
------------------------------------------------------------------------------
$begin set_hes_cross_matrix$$
$spell
$$

$section Setting Hessian Cross Matrix$$

$head Syntax$$
$code%CppAD::mixed::set_hes_cross_matrix(
	%n_fixed%, %n_random%, %row%, %col%, %val%
)%$$

$head Notation$$
We refer to the matrix
$latex f_{u, \theta } ( \theta , u)$$ as the Hessian cross terms.

$head n_fixed$$
This argument has prototype
$codei%
	size_t %n_fixed%
%$$
and is the number of
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$.

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is the number of
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$.
It is assumed that $icode%n_random% > 0%$$.


$subhead row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
It contains one element each possibly non-zero entry in
the Hessian cross terms.
We use the notation $icode%K% = %row%.size()%$$.
For $icode%k% = 0 , %...%, %K%$$,
$icode%row%[%k%]%$$ is the row index of a possibly non-zero entry
in the Hessian cross terms.
Note that for each $icode k$$, $codei%row%[%k%]% < %n_rancom%$$.

$subhead col$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %col%
%$$
and has the same size as $icode row$$.
For $icode%k% = 0 , %...%, %K%$$,
$codei%col%[%k%]%$$ is the column index of a non-zero entry
in the Hessian cross terms.
Note that for each $icode k$$, $codei%col%[%k%]% < %n_fixed%$$.

$subhead val$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val%
%$$
and has the same size as $icode row$$.
For $icode%k% = 0 , %...%, %K%$$,
$icode%val%[%k%]%$$ is the value of a non-zero entry
in the Hessian cross terms.

$subhead hes_cross_matrix_$$
The input value of this matrx does not matter.
Upon return it is an $code Eigen$$ sparse matrix representation
of the Hessian cross terms.

$end
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void set_hes_cross_matrix(
	size_t                       n_fixed  ,
	size_t                       n_random ,
	const CppAD::vector<size_t>& row      ,
	const CppAD::vector<size_t>& col      ,
	const CppAD::vector<size_t>& val      )
{	assert( n_random > 0 );
	assert( row.size() == col.size() );
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
	assert( max_row_p1 <= n_random );
	assert( max_col_p1 <= n_fixed );

	// size of the matrix
	nrow = n_random;
	ncol = n_fixed;
	hes_cross_matrix_.resize(nrow, ncol);

	// put the values in the matrix
	for(size_t k = 0; k < K; k++)
		hes_cross_matrix_.insert(row[k], col[k]) = val[k];

	return;
}
