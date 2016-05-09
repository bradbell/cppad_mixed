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
$begin ldlt_cholmod_sim_cov$$
$spell
	ldlt_obj
	sim_cov
	const
	cholmod
	bool
	xam
	CppAD
$$

$section Simulations with Covariance Corresponding to Factored Matrix$$

$head Syntax$$
$icode%ok% = %ldlt_obj%.sim_cov(%w%, %v%)%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function simulates a normal random vector with mean zero
and covariance $latex H^{-1}$$
where $latex H$$ is the symmetric matrix that has been factored as
$latex \[
	L D L^\R{T} = P H P^\R{T}
\] $$
where $latex L$$ is lower triangular, $latex D$$ is diagonal,
and $latex P$$ is a permutation matrix.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_cholmod %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_cholmod_update$$.

$head w$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %w%
%$$
and its size is equal to the number of rows in $latex H$$.

$head v$$
This argument has prototype
$codei%
	CppAD::vector<double>& %v%
%$$
and its size is equal to the number of rows in $latex H$$.
The input value of its elements does not matter.
Upon return
$latex \[
	v = P^\R{T} L^{-\R{T}} D^{-1/2} w
\] $$
If $latex w$$ is mean zero, variance identity white noise,
$latex w \sim \B{N} ( 0 , I )$$,
then $latex v$$ will be mean zero and variance $latex H^{-1}$$,
$latex v \sim \B{N} ( 0 , H^{-1} )$$; see
$cref/sparse observed information/theory/Sparse Observed Information/$$.

$head ok$$
The return value has prototype
$codei%
	bool %ok%
%$$
It is true if all the elements of $latex D$$ are greater then zero; i.e.,
if $latex H$$ is positive definite.
If $icode ok$$ is false,
$latex H$$ is not positive definite and the values in $latex v$$ are
the same as their input values.

$head Example$$
The file $cref/ldlt_cholmod_xam.cpp/ldlt_cholmod_xam.cpp/sim_cov/$$ contains an
example and test that uses this function.

$end
*/
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/utility/index_sort.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

bool ldlt_cholmod::sim_cov(
	const CppAD::vector<double>& w  ,
	CppAD::vector<double>&       v  )
{	//
	assert( factor_  != CPPAD_NULL );
	assert( rhs_     != CPPAD_NULL );
	assert( rhs_set_ != CPPAD_NULL );
	assert( rhs_->nrow == nrow_ );
	assert( rhs_->ncol == 1     );
	//
	assert( v.size() == nrow_ );
	//
	int*    factor_p  = (int *) factor_->p;
	double* factor_x  = (double *) factor_->x;
	int*    rhs_set_p = (int *) rhs_set_->p;
	int*    rhs_set_i = (int *) rhs_set_->i;
	double* rhs_x     = (double *) rhs_->x;
# ifndef NDEBUG
	int*    factor_i  = (int *) factor_->i;
# endif
	//
	// we will solve for all components on the right hand side
	rhs_set_p[0] = 0;
	rhs_set_p[1] = nrow_;
	for(size_t i = 0; i < nrow_; i++)
		rhs_set_i[i] = i;
	//
	// set rhs_ = w
	for(size_t i = 0; i < nrow_; i++)
		rhs_x[i] = w[i];
	//
	// set rhs_ = D^{-1/2} w
	for(size_t i = 0; i < nrow_; i++)
	{	// first element of each column is always the diagonal element
		assert( size_t( factor_i[ factor_p[i] ] ) == i );
		double di = factor_x[ factor_p[i] ];
		if( di <= 0.0 )
			return false;
		rhs_x[i]   = rhs_x[i] / std::sqrt( di );
	}
	// set sol_ = L^{-T} * D^{-1/2} w
	int sys = CHOLMOD_Lt;
# ifndef NDEBUG
	int flag =
# endif
	cholmod_solve2(
		sys,
		factor_,
		rhs_,
		rhs_set_,
		&sol_,
		&sol_set_,
		&work_one_,
		&work_two_,
		&common_
	);
	// check assumptions
	assert( flag == CHOLMOD_TRUE );
	assert( sol_set_->nrow == nrow_ );
	assert( sol_set_->ncol == 1 );
	assert( sol_set_->xtype == CHOLMOD_PATTERN );
	assert( sol_set_->packed == CHOLMOD_TRUE);
# ifndef NDEBUG
	int*   sol_set_p = (int *) sol_set_->p;
	assert( size_t( sol_set_p[1] ) == nrow_ );
# endif
	assert( sol_ != CPPAD_NULL   );
	assert( sol_->nrow == nrow_ );
	assert( sol_->ncol == 1     );
	//
	// set rhs_ = L^{-T} * D^{-1/2} w
	double* sol_x = (double *) sol_->x;
	for(size_t i = 0; i < nrow_; i++)
		rhs_x[i] = sol_x[i];
	//
	// set sol_ = P^T L^{-T} * D^{-1/2} w
	sys = CHOLMOD_Pt;
# ifndef NDEBUG
	flag =
# endif
	cholmod_solve2(
		sys,
		factor_,
		rhs_,
		rhs_set_,
		&sol_,
		&sol_set_,
		&work_one_,
		&work_two_,
		&common_
	);
	// check assumptions
	assert( flag == CHOLMOD_TRUE );
	assert( sol_set_->nrow == nrow_ );
	assert( sol_set_->ncol == 1 );
	assert( sol_set_->xtype == CHOLMOD_PATTERN );
	assert( sol_set_->packed == CHOLMOD_TRUE);
# ifndef NDEBUG
	sol_set_p = (int *) sol_set_->p;
	assert( size_t( sol_set_p[1] ) == nrow_ );
# endif
	assert( sol_ != CPPAD_NULL   );
	assert( sol_->nrow == nrow_ );
	assert( sol_->ncol == 1     );
	//
	// return v
	sol_x = (double *) sol_->x;
	for(size_t i = 0 ; i < nrow_; i++)
		v[i] = sol_x[i];
	//
	return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
