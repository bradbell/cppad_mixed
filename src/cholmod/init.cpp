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
------------------------------------------------------------------------------
$begin cholmod_init$$
$spell
	Cholesky
	chol_ran_hes
	CppAD
	const
	init
	cholmod
	pos
	rhs
$$

$section Initialize Cholesky Factor for a Specific Sparsity Pattern$$

$head Syntax$$
$icode%chol_ran_hes%.init(%n_fixed%, %n_random%, %hes_info%)
%$$

$head Private$$
The $cref cholmod$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head chol_ran_hes$$
This object has prototype
$codei%
	CppAD::mixed::cholmod %chol_ran_hes%
%$$

$head n_fixed$$
This argument has prototype
$codei%
	size_t %n_fixed%
%$$
and is then number of fixed effects.

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is then number of random effects.

$head hes_info$$
This argument had prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It is a
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
for the Hessian of the
$cref/random likelihood/theory/Random Likelihood, f(theta, u)/$$
$latex f_{u,u} ( \theta , u )$$.
Furthermore, it is in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order.

$subhead row$$
For $icode%k% = 0 , %...% , %hes_info%.row.size()-1%$$,
$codei%
	%n_fixed% <= %row%[%k%] < %n_fixed% + %n_random%
%$$
The reason for the offset is that these indices are relative to both
the fixed and random effects and the fixed effects come before the
random effects.

$subhead col$$
For $icode%k% = 0 , %...% , %hes_info%.col%.size()-1%$$,
$codei%
	%n_fixed% <= %col%[%k%] <= %row%[%k%]
%$$

$head Assumptions$$
All of the $code cholmod$$ private pointers
are null when this routine is called.

$head pos_matrix_$$
This is set to a packed real lower triangular real sparse matrix
with the pattern specified by $icode hes_info$$ and
the value $code nan$$ for each possibly non-zero value.

$head factor_$$
This is the result of
$codei%
	factor_ = cholmod_analyze(pos_matrix_, &common_)
%$$

$head rhs$$
This is set to dense column vector of $code n_random$$ zeros as follows:
$codei%
	rhs_ = cholmod_zeros(%n_random%, 1, CHOLMOD_REAL, &common_)
%$$

$end
*/

# include <cppad/mixed/cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cholmod::init(
	size_t                               n_fixed  ,
	size_t                               n_random ,
	const CppAD::mixed::sparse_mat_info& hes_info )
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
	size_t nrow  = n_random;
	size_t ncol  = n_random;
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
	{	assert( n_fixed <= hes_info.row[k] );
		assert( n_fixed <= hes_info.col[k] );
		T_i[k] = static_cast<int>( hes_info.row[k] - n_fixed );
		T_j[k] = static_cast<int>( hes_info.col[k] - n_fixed );
		T_x[k] = nan;
	}
	triplet->nnz = nzmax;

	// convert triplet to a sparse matrix
	nzmax = 0; // just need max corresponding to the triplet
	pos_matrix_ = cholmod_triplet_to_sparse(triplet, nzmax, &common_);

	// check assumptions
	assert( pos_matrix_->nrow   == n_random );
	assert( pos_matrix_->ncol   == n_random );
	assert( pos_matrix_->nzmax  == hes_info.row.size() );
	assert( pos_matrix_->stype  == CHOLMOD_STYPE_LOWER_TRIANGLE );
	assert( pos_matrix_->itype  == CHOLMOD_INT );
	assert( pos_matrix_->xtype  == CHOLMOD_REAL );
	assert( pos_matrix_->dtype  == CHOLMOD_DOUBLE );
	assert( pos_matrix_->sorted == CHOLMOD_FALSE );
	assert( pos_matrix_->packed == CHOLMOD_TRUE  );


	// analyze the sparsity pattern for LDL^T Cholesky factorization of
	// f_{u,u}(theta, u)
	factor_ = cholmod_analyze(pos_matrix_, &common_);

	// set rhs_ to column vector of zeros
	rhs_ = cholmod_zeros(n_random, 1, CHOLMOD_REAL, &common_);

	// done with triplet
	cholmod_free_triplet(&triplet, &common_ );
}

} } // END_CPPAD_MIXED_NAMESPACE
