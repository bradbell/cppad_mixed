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
$begin cholmod_init$$
$spell
	nrow
	xam
	Cholesky
	CppAD
	const
	init
	cholmod_obj
	pos
	rhs
	hes
$$

$section Initialize Cholesky Factor for a Specific Sparsity Pattern$$

$head Syntax$$
$icode%cholmod_obj%.init(%hes_info%)
%$$

$head Private$$
The $cref cholmod$$ class is an
$cref/implementation detail/choleig/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head cholmod_obj$$
This object has prototype
$codei%
	CppAD::mixed::cholmod %cholmod_obj%
%$$

$head hes_info$$
This argument had prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It is a
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$ for the
square matrices with
$cref/nrow_/cholmod_ctor/nrow_/$$ rows that we will compute the
Cholesky factor of.
It is in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$head Assumptions$$
All of the $code cholmod$$ private pointers
are null when this routine is called.

$head pos_matrix_$$
Upon return,
this is set to a packed, real, sorted, lower triangular, sparse matrix
with the pattern specified by $icode hes_info$$ and
the value $code nan$$ for each possibly non-zero value.

$head factor_$$
Upon return,
this is the result of
$codei%
	factor_ = cholmod_analyze(pos_matrix_, &common_)
%$$

$head rhs$$
Upon return,
this is set to dense column vector of $code nrow_$$ zeros as follows:
$codei%
	rhs_ = cholmod_zeros(nrow_, 1, CHOLMOD_REAL, &common_)
%$$

$head Example$$
The file $cref/cholmod_xam.cpp/cholmod_xam.cpp/init/$$ contains an
example and test that uses this function.

$end
*/

# include <cppad/mixed/cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cholmod::init( const CppAD::mixed::sparse_mat_info& hes_info )
{
	assert(pos_matrix_ == CPPAD_NULL );
	assert(factor_     == CPPAD_NULL );
	assert(rhs_        == CPPAD_NULL );
	assert(rhs_set_    == CPPAD_NULL );
	assert(sol_        == CPPAD_NULL );
	assert(sol_set_    == CPPAD_NULL );
	assert(work_one_   == CPPAD_NULL );
	assert(work_two_   == CPPAD_NULL );

	// set triplet corresponding to sparsity pattern and reserve spase
	// for unspecified values (that end up in pos_matrix_.
	size_t nrow  = nrow_;
	size_t ncol  = nrow_;
	size_t nzmax = hes_info.row.size();
	int    stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
	int    xtype = CHOLMOD_REAL;
	cholmod_triplet* triplet = cholmod_allocate_triplet(
		nrow, ncol, nzmax, stype, xtype, &common_
	);
	int* T_i    = (int *) triplet->i;
	int* T_j    = (int *) triplet->j;
	double* T_x = (double *) triplet->x;
	double nan  = std::numeric_limits<double>::quiet_NaN();
	for(size_t k = 0; k < nzmax; k++)
	{	assert( hes_info.row[k] < nrow_ );
		assert( hes_info.col[k] < nrow_ );
		//
		T_i[k] = static_cast<int>( hes_info.row[k] );
		T_j[k] = static_cast<int>( hes_info.col[k] );
		//
		T_x[k] = nan;
	}
	triplet->nnz = nzmax;

	// convert triplet to a sparse matrix
	nzmax = 0; // just need max corresponding to the triplet
	pos_matrix_ = cholmod_triplet_to_sparse(triplet, nzmax, &common_);

	// check assumptions
	assert( pos_matrix_->nrow   == nrow_ );
	assert( pos_matrix_->ncol   == nrow_ );
	assert( pos_matrix_->nzmax  == hes_info.row.size() );
	assert( pos_matrix_->stype  == CHOLMOD_STYPE_LOWER_TRIANGLE );
	assert( pos_matrix_->itype  == CHOLMOD_INT );
	assert( pos_matrix_->xtype  == CHOLMOD_REAL );
	assert( pos_matrix_->dtype  == CHOLMOD_DOUBLE );
	assert( pos_matrix_->sorted == CHOLMOD_TRUE  );
	assert( pos_matrix_->packed == CHOLMOD_TRUE  );


	// analyze the sparsity pattern for LDL^T Cholesky factorization of
	factor_ = cholmod_analyze(pos_matrix_, &common_);

	// check assumptions
	assert( factor_->n     == nrow_ );
	assert( factor_->minor == nrow_ );
	assert( factor_->is_ll == CHOLMOD_FALSE );

	// set rhs_ to column vector of zeros
	rhs_ = cholmod_zeros(nrow_, 1, CHOLMOD_REAL, &common_);

	// done with triplet
	cholmod_free_triplet(&triplet, &common_ );
}

} } // END_CPPAD_MIXED_NAMESPACE
