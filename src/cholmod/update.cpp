/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ldlt_cholmod_update$$
$spell
	rcv
	bool
	sym
	ldlt
	xam
	const
	CppAD
	hes
	ldlt_obj
	cholmod
	init
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section Update Factorization Using new Matrix Values$$

$head Syntax$$
$icode%ok% = %ldlt_obj%.update(%H_rcv%)%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
CppAD Mixed user API.

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

$head H_rcv$$
This argument contains new values for the
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
we are computing the LDLT factor of.
The $cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
must be the same as in $cref/ldlt_cholmod_init/ldlt_cholmod_init/H_rc/$$.
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
using the sparsity pattern as is in $icode H_rcv$$.
Upon return, values in $code sym_matrix_$$ have been set to the
corresponding values in the vector $icode%H_rcv%.val()%$$.

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

$head ok$$
If $icode ok$$ is true, the matrix was factored.
Otherwise, the matrix is singular.

$head Order of Operations$$
This $icode ldlt_obj$$ function must be called,
after the constructor and $cref/init/ldlt_cholmod_init/$$
and before any other member functions
(except $cref/pattern/ldlt_cholmod_pattern/$$).

$head Example$$
The file $cref/ldlt_cholmod.cpp/ldlt_cholmod.cpp/update/$$ contains an
example and test that uses this function.

$end
*/
// ----------------------------------------------------------------------------


# include <cppad/mixed/include_cholmod.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// BEGIN_PROTOTYPE
bool ldlt_cholmod::update(const CppAD::mixed::d_sparse_rcv& H_rcv)
// END_PROTOTYPE
{	assert( init_done_ );
	update_called_ = true;
	//
	// set the values in sym_matrix_
	double* H_x  = reinterpret_cast<double*>( sym_matrix_->x);
# ifndef NDEBUG
	int*    H_p  = reinterpret_cast<int *>( sym_matrix_->p );
	int*    H_i  = reinterpret_cast<int *>( sym_matrix_->i );
# endif
	for(size_t k = 0; k < H_rcv.nnz(); k++)
	{	size_t ell = H_rc2cholmod_order_[k];
		H_x[ell]   = H_rcv.val()[k];
# ifndef NDEBUG
		assert( size_t( H_i[ell] ) == H_rcv.row()[k] );
		size_t j = H_rcv.col()[k];
		assert( size_t(H_p[j]) <= ell && ell < size_t(H_p[j+1]) );
# endif
	}
	// set factor_ to LDL^T factorization for this value of Hessian
# ifndef NDEBUG
	int flag =
# endif
	cholmod_factorize(sym_matrix_, factor_, &common_);
	//
	if( common_.status == CHOLMOD_NOT_POSDEF )
		return false;
	else
		assert( common_.status == CHOLMOD_OK );

	// check assumptions
	assert( flag           == CHOLMOD_TRUE );
	assert( factor_->n     == nrow_ );
	assert( factor_->minor <= nrow_ );
	assert( factor_->is_ll == CHOLMOD_FALSE );
	assert( factor_->xtype == CHOLMOD_REAL );
	// ---------------------------------------------------------------------
	if( factor_->minor < nrow_ )
		return false;
	return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
