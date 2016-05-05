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
$begin ldlt_cholmod_update$$
$spell
	sym
	ldlt
	Cholesky
	xam
	const
	CppAD
	hes
	ldlt_obj
	cholmod
	init
	pos
$$

$section Update LDLT Factorization Using new Matrix Values$$

$head Syntax$$
$icode%pos% = ldlt_obj%.update(%hes_info%)%$$


$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This routine updates the $cref ldlt_cholmod$$ factorization
for new values in the square symmetric matrix.

$head ldlt_obj$$
This object has prototype
$codei%
	CppAD::mixed::ldlt_cholmod %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_cholmod_init$$.

$head hes_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It contains new values for the
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
we are computing the Cholesky factor of.
The $cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
must be the same as in $cref/ldlt_cholmod_init/ldlt_cholmod_init/hes_info/$$.
Hence, in particular, it must be in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$head sym_matrix_$$
On input, the member variable
$codei%
	cholmod_sparse* sym_matrix_
%$$
has been $cref/initialized/ldlt_cholmod_init/sym_matrix_/$$
using the sparsity pattern as is in $icode hes_info$$.
Upon return, values in $code sym_matrix_$$ have been set to the
corresponding values in the vector $icode%hes_info%.val%$$.

$head factor_$$
On input, the member variable
$codei%
	cholmod_factor* factor_
%$$
has been $cref/initialized/ldlt_cholmod_init/factor_/$$
using the sparsity pattern for the Hessian.
Upon return, it contains the factorization
$codei%
	cholmod_factorize(sym_matrix_, factor_, &common_)
%$$

$head pos$$
If the matrix is positive definite, $icode pos$$ is true.
Otherwise, it is false.

$head Example$$
The file $cref/ldlt_cholmod_xam.cpp/ldlt_cholmod_xam.cpp/update/$$ contains an
example and test that uses this function.

$end
*/
// ----------------------------------------------------------------------------


# include <cholmod.h>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

bool ldlt_cholmod::update( const CppAD::mixed::sparse_mat_info& hes_info )
{
	// set the values in sym_matrix_
	double* H_x = (double *) sym_matrix_->x;
	int*    H_p = (int *)    sym_matrix_->p;
# ifndef NDEBUG
	int*    H_i = (int *)    sym_matrix_->i;
# endif
	for(size_t j = 0; j < nrow_; j++)
	{	for(size_t k = (size_t) H_p[j]; k < (size_t) H_p[j+1]; k++)
		{	assert( hes_info.row[k] == (size_t) H_i[k] );
			assert( hes_info.col[k] == j );
			H_x[k] = hes_info.val[k];
		}
	}
	// set factor_ to LDL^T factorization for this value of Hessian
# ifndef NDEBUG
	int flag =
# endif
	cholmod_factorize(sym_matrix_, factor_, &common_);

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
	return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
