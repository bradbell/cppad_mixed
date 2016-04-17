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
$begin cholmod_update$$
$spell
	Cholesky
	xam
	const
	CppAD
	hes
	cholmod_obj
	init
	pos
$$

$section Update Factorization Using new Matrix Values$$

$head Syntax$$
$icode%pos% = cholmod_obj%.update(%hes_info%)%$$


$head Private$$
The $cref cholmod$$ class is an
$cref/implementation detail/cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This routine updates the $cref cholmod$$ factorization
for new values in the square positive definite matrix.

$head cholmod_obj$$
This object has prototype
$codei%
	CppAD::mixed::cholmod %cholmod_obj%
%$$
In addition, it must have a previous call to
$cref cholmod_init$$.

$head hes_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It contains new values for the
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
we are computing the Cholesky factor of.
The $cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
must be the same as in $cref/cholmod_init/cholmod_init/hes_info/$$.
Hence, in particular, it must be in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$head pos_matrix_$$
On input, the member variable
$codei%
	cholmod_sparse* pos_matrix_
%$$
has been $cref/initialized/cholmod_init/pos_matrix_/$$
using the sparsity pattern as is in $icode hes_info$$.
Upon return, values in $code pos_matrix_$$ have been set to the
corresponding values in the vector $icode%hes_info%.val%$$.

$head factor_$$
On input, the member variable
$codei%
	cholmod_factor* factor_
%$$
has been $cref/initialized/cholmod_init/factor_/$$
using the sparsity pattern for the Hessian.
Upon return, it contains the factorization
$codei%
	cholmod_factorize(pos_matrix_, factor_, &common_)
%$$

$head pos$$
If the matrix is positive definite, $icode pos$$ is true.
Otherwise, it is false.

$head Example$$
The file $cref/cholmod_xam.cpp/cholmod_xam.cpp/update/$$ contains an
example and test that uses this function.

$end
*/
// ----------------------------------------------------------------------------


# include <cholmod.h>
# include <cppad/mixed/cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

bool cholmod::update( const CppAD::mixed::sparse_mat_info& hes_info )
{
	// set the values in pos_matrix_
	double* A_x = (double *) pos_matrix_->x;
	int*    A_p = (int *)    pos_matrix_->p;
# ifndef NDEBUG
	int*    A_i = (int *)    pos_matrix_->i;
# endif
	for(size_t j = 0; j < nrow_; j++)
	{	for(size_t k = (size_t) A_p[j]; k < (size_t) A_p[j+1]; k++)
		{	assert( hes_info.row[k] == (size_t) A_i[k] );
			assert( hes_info.col[k] == j );
			A_x[k] = hes_info.val[k];
		}
	}
	// set factor_ to LDL^T factorization for this value of Hessian
# ifndef NDEBUG
	int flag =
# endif
	cholmod_factorize(pos_matrix_, factor_, &common_);

	// check assumptions
	assert( flag           == CHOLMOD_TRUE );
	assert( factor_->n     == nrow_ );
	assert( factor_->minor == nrow_ );
	assert( factor_->is_ll == CHOLMOD_FALSE );
	assert( factor_->xtype == CHOLMOD_REAL );
	//
	if( common_.status == CHOLMOD_NOT_POSDEF )
		return false;
	assert( common_.status == CHOLMOD_OK );
	// ---------------------------------------------------------------------
	// cholmod_solve2 misses some cases, see test_more/cholmod_solve2.cpp
	int*    L_p  = (int *) factor_->p;
	double* L_x  = (double *) factor_->x;
# ifndef NDEBUG
	int*    L_i  = (int *) factor_->i;
# endif
	for(size_t j = 0; j < nrow_; j++)
	{	// first element for each column is always the diagonal element
		assert( size_t( L_i [ L_p[j] ] ) == j );
		// j-th element on diagonal of D in factorization
		if( L_x[ L_p[j] ] <= 0.0 )
			return false;
	}
	// ---------------------------------------------------------------------
	return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
