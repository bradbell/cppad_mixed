// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cholmod.h>
# include <cppad/mixed/cholmod.hpp>
/*
$begin cholmod_update$$
$spell
	const
	CppAD
	hes
	cholmod
	init
	pos
	chol
$$

$section Update the Factorization of Hessian w.r.t. Random Effects$$

$head Syntax$$
$icode%chol_ran_hes%.cholmod_update(%hes_info%)%$$


$head Private$$
The $cref cholmod$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This routine updates the $cref cholmod$$ factorization
for new values in the Hessian  $latex f_{u,u} ( \theta , u )$$.

$head chol_ran_hes$$
This is a $cref cholmod$$ object.

$head hes_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It contains a new
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$ representation
of the Hessian $latex f_{u,u} ( \theta , u )$$.
The $cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
is must be the same as in $cref/cholmod_init/cholmod_init/hes_info/$$.
Hence, in particular, it must be in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order.
In addition, the row and column indices in $code pos_matrix_$$ are
shifted by $code n_fixed_$$.

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

$comment%
	example/private/cholmod_update.cpp
%$$
$head Example$$
The file $code cholmod_update.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------


# include <cppad/mixed/cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cholmod::update( const CppAD::mixed::sparse_mat_info& hes_info )
{
	// set the values in pos_matrix_
	double* A_x = (double *) pos_matrix_->x;
	int*    A_p = (int *)    pos_matrix_->p;
	int*    A_i = (int *)    pos_matrix_->i;
	for(size_t j = 0; j < n_random_; j++)
	{	for(size_t k = (size_t) A_p[j]; k < (size_t) A_p[j+1]; k++)
		{	assert( hes_info.row[k] == n_fixed_ + (size_t) A_i[k] );
			assert( hes_info.col[k] == n_fixed_ + j );
			A_x[k] = hes_info.val[k];
		}
	}
	// set factor_ to LDL^T factorization for this value of Hessian
	cholmod_factorize(pos_matrix_, factor_, &common_);
}

} } // END_CPPAD_MIXED_NAMESPACE
