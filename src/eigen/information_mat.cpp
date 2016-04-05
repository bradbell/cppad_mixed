// $Id:$
/*
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
$begin information_mat$$
$spell
	CppAD
	cppad
$$

$section Compute the Observed Information Matrix$$

$head Syntax$$
$icode%information_info% = %mixed_object%.information_mat(
	%solution%, %random_options%, %random_lower%, %random_upper%, %random_in%
)%$$

$head Purpose$$
Compute the observed information matrix.
We use $latex L ( \theta )$$ to denote the
$cref/total objective/theory/Objective/Total Objective, L(theta)/$$.
If $latex \hat{\theta}$$ is the optimal value for the fixed effects,
the corresponding observed information is
$latex \[
	L^{(2)} ( \hat{\theta} )
\]$$
Absolute value terms in the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
for the $cref fix_likelihood$$ are not include in this Hessian
(because they do not have a derivative, let alone Hessian, at zero).

$head quasi_fixed$$
If $cref/quasi_fixed/derived_ctor/quasi_fixed/$$ was true when
$icode mixed_object$$ was constructed,
The $cref initialize$$ routine did not include the full Newton method
(hence uses less memory).
This memory is allocated and used during $code sample_fixed$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head solution$$
is the $cref/solution/optimize_fixed/solution/$$
for a previous call to $cref optimize_fixed$$.

$head random_options$$
is the $cref/random_options/optimize_fixed/random_options/$$
for a previous call to $code optimize_fixed$$.

$head random_lower$$
is the $cref/random_lower/optimize_fixed/random_lower/$$
for a previous call to $code optimize_fixed$$.

$head random_upper$$
is the $cref/random_upper/optimize_fixed/random_upper/$$
for a previous call to $code optimize_fixed$$.

$head random_in$$
is the $cref/random_in/optimize_fixed/random_in/$$
for a previous call to $code optimize_fixed$$.

$head information_info$$
The return value has prototype
$codei%
	CppAD::mixed::sparse_mat_info %information_info%
%$$
This is a sparse matrix representation for the
lower triangle of the observed information matrix.
(Note that this matrix is symmetric and hence determined by its
lower triangle.)

$end
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>

CppAD::mixed::sparse_mat_info cppad_mixed::information_mat(
	const CppAD::mixed::fixed_solution&  solution             ,
	const std::string&                   random_options       ,
	const d_vector&                      random_lower         ,
	const d_vector&                      random_upper         ,
	const d_vector&                      random_in            )
{	using Eigen::Dynamic;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
	typedef eigen_sparse::InnerIterator                       sparse_itr;
	// solution
	assert( solution.fixed_opt.size()   == n_fixed_ );
	assert( solution.fixed_lag.size()   == n_fixed_ );
	// random_(lower, upper, in)
	assert( random_lower.size() == n_random_ );
	assert( random_upper.size() == n_random_ );
	assert( random_in.size()    == n_random_ );
	//
	// optimal fixed effects
	const d_vector& fixed_opt( solution.fixed_opt );
	// -----------------------------------------------------------------------
	// optimize the random effects
	d_vector random_opt = optimize_random(
		random_options, fixed_opt, random_lower, random_upper, random_in
	);
	// update the cholesky factor for this fixed and random effect
	update_factor(fixed_opt, random_opt);
	// -----------------------------------------------------------------------
	// If Quasi-Newton method was used, must initilaize routines
	// that are only used for the Hessian calculation; see initilaize.cpp
	if( n_random_ != 0 && ! init_newton_atom_done_ )
	{	assert( quasi_fixed_ );
		assert( ! init_ran_objcon_done_ );
		assert( ! init_ran_objcon_hes_done_ );
		//
		// newton_atom_
		assert( ran_like_a1fun_.size_var() > 0  );
		newton_atom_.initialize(
			ran_like_a1fun_, fixed_opt, random_opt
		);
		assert( init_newton_atom_done_ );
		//
		// ran_objcon_fun_
		assert( ! init_ran_objcon_done_ );
		init_ran_objcon(fixed_opt, random_opt);
		assert( init_ran_objcon_done_ );
		//
		// ran_objcon_hes_
		assert( ! init_ran_objcon_hes_done_ );
		init_ran_objcon_hes(fixed_opt, random_opt);
		assert( init_ran_objcon_hes_done_ );
	}
	assert( init_newton_atom_done_ );
	assert( init_ran_objcon_done_ );
	assert( init_ran_objcon_hes_done_ );
	// -----------------------------------------------------------------------
	// Lower triangle of Hessian w.r.t. fixed effects
	// for random part of objective,no constraints
	d_vector w_ran(n_ran_con_ + 1);
	w_ran[0] = 1.0;
	for(size_t j = 1; j <=n_ran_con_; j++)
		w_ran[j] = 0.0;
	//
	CppAD::mixed::sparse_mat_info ran_info;
	ran_objcon_hes(
		fixed_opt, random_opt, w_ran, ran_info.row, ran_info.col, ran_info.val
	);
	eigen_sparse ran_hes = CppAD::mixed::triple2eigen(
		n_fixed_      ,
		n_fixed_      ,
		ran_info.row  ,
		ran_info.col  ,
		ran_info.val
	);
	// -----------------------------------------------------------------------
	// Lower triangle of Hessian of the fixed likelihood
	size_t n_fix_like = 0;
	if( fix_like_fun_.size_var() != 0 )
		n_fix_like = fix_like_fun_.Range();
	d_vector w_fix(n_fix_like);
	w_fix[0] = 1.0;
	for(size_t j = 1; j < n_fix_like; j++)
		w_fix[0] = 0.0;
	//
	CppAD::mixed::sparse_mat_info fix_info;
	fix_like_hes(
			fixed_opt, w_fix, fix_info.row, fix_info.col, fix_info.val
	);
	eigen_sparse fix_hes = CppAD::mixed::triple2eigen(
		n_fixed_      ,
		n_fixed_      ,
		fix_info.row  ,
		fix_info.col  ,
		fix_info.val
	);
	// -----------------------------------------------------------------------
	// Hessian of total objective (observed information matrix)
	eigen_sparse total_hes = ran_hes + fix_hes;
	CppAD::mixed::sparse_mat_info total_info;
	for(size_t j = 0; j < n_fixed_; j++)
	{	for(sparse_itr itr(total_hes, j); itr; ++itr)
		{	assert( size_t( itr.col() ) == j );
			total_info.row.push_back( itr.row() );
			total_info.col.push_back( itr.col() );
			total_info.val.push_back( itr.value() );
		}
	}
	return total_info;
}