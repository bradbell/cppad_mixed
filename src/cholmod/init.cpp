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
$begin ldlt_cholmod_init$$
$spell
	sym
	ldlt
	nrow
	xam
	CppAD
	const
	init
	ldlt_obj
	cholmod
	pos
	rhs
	hes
$$

$section Initialize Factor for a Specific Sparsity Pattern$$

$head Syntax$$
$icode%ldlt_obj%.init(%H_info%)
%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head ldlt_obj$$
This object has prototype
$codei%
	CppAD::mixed::ldlt_cholmod %ldlt_obj%
%$$

$head H_info$$
This argument had prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %H_info%
%$$
It is a
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$ for the
square matrices with
$cref/nrow_/ldlt_cholmod_ctor/nrow_/$$ rows that we will compute the
LDLT factor of.
It is in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$head Assumptions$$
All of the $code cholmod$$ private pointers
are null when this routine is called.

$head sym_matrix_$$
Upon return,
this is set to a packed, real, sorted, lower triangular, sparse matrix
with the pattern specified by $icode H_info$$ and
the value $code nan$$ for each possibly non-zero value.

$head factor_$$
Upon return,
this is the result of
$codei%
	factor_ = cholmod_analyze(sym_matrix_, &common_)
%$$

$head rhs_$$
Upon return,
this is set to dense column vector of $code nrow_$$ zeros as follows:
$codei%
	rhs_ = cholmod_zeros(nrow_, 1, CHOLMOD_REAL, &common_)
%$$

$head rhs_set_$$
Upon return,
this is a sparse column vector with $code nrow_$$ rows.
It is packed, sorted, not symmetric, and just a sparsity pattern.
(Note that sorted refers the column order and does not matter for vectors.)
There are $code nrow_$$
possibly non-zero elements, but which elements
are non-zero is not specified
(expected to be set for each $cref ldlt_cholmod_solve_H$$ usage).

$head Example$$
The file $cref/ldlt_cholmod_xam.cpp/ldlt_cholmod_xam.cpp/init/$$ contains an
example and test that uses this function.

$end
*/

# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/utility/index_sort.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void ldlt_cholmod::init( const CppAD::mixed::sparse_mat_info& H_info )
{
	assert(sym_matrix_ == CPPAD_NULL );
	assert(factor_     == CPPAD_NULL );
	assert(rhs_        == CPPAD_NULL );
	assert(rhs_set_    == CPPAD_NULL );
	assert(sol_        == CPPAD_NULL );
	assert(sol_set_    == CPPAD_NULL );
	assert(work_one_   == CPPAD_NULL );
	assert(work_two_   == CPPAD_NULL );

	// set triplet corresponding to sparsity pattern and reserve spase
	// for unspecified values (that end up in sym_matrix_.
	size_t ncol  = nrow_;
	size_t nzmax = H_info.row.size();
	int    stype = CHOLMOD_STYPE_LOWER_TRIANGLE;
	int    xtype = CHOLMOD_REAL;
	cholmod_triplet* triplet = cholmod_allocate_triplet(
		nrow_, ncol, nzmax, stype, xtype, &common_
	);
	int* T_i    = (int *) triplet->i;
	int* T_j    = (int *) triplet->j;
	double* T_x = (double *) triplet->x;
	double nan  = std::numeric_limits<double>::quiet_NaN();
	for(size_t k = 0; k < nzmax; k++)
	{	assert( H_info.row[k] < nrow_ );
		assert( H_info.col[k] < nrow_ );
		//
		T_i[k] = static_cast<int>( H_info.row[k] );
		T_j[k] = static_cast<int>( H_info.col[k] );
		//
		T_x[k] = nan;
	}
	triplet->nnz = nzmax;

	// convert triplet to a sparse matrix
	nzmax = 0; // just need max corresponding to the triplet
	sym_matrix_ = cholmod_triplet_to_sparse(triplet, nzmax, &common_);

	// check assumptions
	assert( sym_matrix_->nrow   == nrow_ );
	assert( sym_matrix_->ncol   == nrow_ );
	assert( sym_matrix_->nzmax  == H_info.row.size() );
	assert( sym_matrix_->stype  == CHOLMOD_STYPE_LOWER_TRIANGLE );
	assert( sym_matrix_->itype  == CHOLMOD_INT );
	assert( sym_matrix_->xtype  == CHOLMOD_REAL );
	assert( sym_matrix_->dtype  == CHOLMOD_DOUBLE );
	assert( sym_matrix_->sorted == CHOLMOD_TRUE  );
	assert( sym_matrix_->packed == CHOLMOD_TRUE  );


	// analyze the sparsity pattern for LDLT factorization of
	factor_ = cholmod_analyze(sym_matrix_, &common_);

	// check assumptions
	assert( factor_->n     == nrow_ );
	assert( factor_->minor == nrow_ );
	assert( factor_->is_ll == CHOLMOD_FALSE );

	// set rhs_ to column vector of zeros
	rhs_ = cholmod_zeros(nrow_, 1, CHOLMOD_REAL, &common_);

	// set rhs_set_ to be a sparsity pattern, but do not specify the
	// actual non-zero entries.
	ncol       = 1;
	nzmax      = nrow_;
	int sorted = CHOLMOD_TRUE;
	int packed = CHOLMOD_TRUE;
	stype      = CHOLMOD_STYPE_NOT_SYMMETRIC;
	xtype      = CHOLMOD_PATTERN;
	rhs_set_   = cholmod_allocate_sparse(
		nrow_,
		ncol,
		nzmax,
		sorted,
		packed,
		stype,
		xtype,
		&common_
	);
	// check
	assert( rhs_set_->nrow   == nrow_ );
	assert( rhs_set_->ncol   == 1     );
	assert( rhs_set_->nzmax  == nrow_ );
	assert( rhs_set_->stype  == CHOLMOD_STYPE_NOT_SYMMETRIC );
	assert( rhs_set_->itype  == CHOLMOD_INT );
	assert( rhs_set_->sorted == CHOLMOD_TRUE  );
	assert( rhs_set_->packed == CHOLMOD_TRUE  );

	// done with triplet
	cholmod_free_triplet(&triplet, &common_ );

	// set info2cholmod_order_
	int*    H_p = (int *)    sym_matrix_->p;
	int*    H_i = (int *)    sym_matrix_->i;
	assert( size_t( H_p[nrow_] ) == H_info.row.size() );
	nzmax = H_info.row.size();
	CppAD::vector<size_t> key(nzmax);
	info2cholmod_order_.resize(nzmax);
	size_t k = 0;
	ncol     = nrow_;
	for(size_t j = 0; j < ncol; j++)
	{	size_t start = size_t( H_p[j] );
		size_t end   = size_t( H_p[j+1] );
		for(size_t ell = start; ell < end; ell++)
		{	// H_info is in column major order
			assert( H_info.col[k] == j );
			size_t i = size_t(H_i[ell]);
			assert( i < nrow_ );
			// column major order for sym_martrix_
			key[k++] = j * nrow_ + i;
		}
	}
	assert( k == nzmax );
	CppAD::index_sort(key, info2cholmod_order_);
# ifndef NDEBUG
	for(size_t k = 0; k < nzmax; k++)
		assert( size_t(H_i[ info2cholmod_order_[k] ]) == H_info.row[k] );
# endif
}

} } // END_CPPAD_MIXED_NAMESPACE
