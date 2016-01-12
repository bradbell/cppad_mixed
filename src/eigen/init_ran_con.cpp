// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin init_ran_con$$
$spell
	init
	cppad
	const
	CppAD
	cholesky
	eigen
$$

$section Initialize Random Constraints$$

$head Syntax$$
$icode%mixed_object%.init_ran_con(%n_random%, %A_info%)
%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is the number of
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$.
It is assumed that $icode%n_random% > 0%$$.

$head A_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %A_info%
%$$
It is a
$cref/sparse matrix/sparse_mat_info/val/Sparse Matrix/$$ representation
of the
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$.

$subhead Restriction$$
We assume that $icode%A_info.val[%k%] != 0%$$ for
$icode%k% = 0 , ... , %K%-1%$$ where
$icode%K% = %A_info%.row.size()%$$.

$head ran_con_mat_$$
The input value of the member variable
$codei%
	CppAD::mixed::cholesky::eigen_sparse ran_con_mat_
%$$
does not matter.
Upon return it is an $code Eigen$$ sparse matrix representation
of the random constraint matrix $latex A$$.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_ran_con(
	size_t                               n_random ,
	const CppAD::mixed::sparse_mat_info& A_info   )
{	assert( ! init_ran_con_done_ );

	// number of possibly non-zero entries
	size_t K = A_info.row.size();
	//
	assert( n_random > 0 );
	assert( A_info.col.size() == K );
	assert( A_info.val.size() == K );

	// determine maximum row and column index
	size_t max_row_p1 = 0;
	size_t max_col_p1 = 0;
	for(size_t k = 0; k < K; k++)
	{	assert( A_info.val[k] != 0.0 );
		max_row_p1 = std::max(max_row_p1, A_info.row[k] + 1);
		max_col_p1 = std::max(max_col_p1,    A_info.col[k] + 1);
	}
	assert( max_col_p1 <= n_random );

	// size of the matrix A
	size_t nrow = max_row_p1;
	size_t ncol = n_random;
	ran_con_mat_.resize(nrow, ncol);
	// we know number of non-zeros
	ran_con_mat_.reserve(K);

	// put the values in A
	for(size_t k = 0; k < K; k++)
		ran_con_mat_.insert(A_info.row[k], A_info.col[k]) = A_info.val[k];
	//
	init_ran_con_done_ = true;
	return;
}
