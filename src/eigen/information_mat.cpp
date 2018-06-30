// $Id:$
/*
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-----------------------------------------------------------------------------
$begin information_mat$$
$spell
	CppAD
	cppad
	rcv
$$

$section Compute the Observed Information For Fixed Effects$$

$head Syntax$$
$icode%information_rcv% = %mixed_object%.information_mat(
	%solution%, %random_opt%
)%$$
$icode%information_rcv% = %mixed_object%.information_mat(
	%solution%, %random_opt%
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

$head random_opt$$
is the optimal random effects corresponding to the solution; i.e.
$codei%
	%random_opt% = %mixed_object%.optimize_random(
		%random_options%,
		%solution%.fixed_opt,
		%random_lower%,
		%random_upper%,
		%random_in%
	)
%$$
$icode random_options$$,
$icode random_lower$$,
$icode random_upper$$, and
$icode random_in$$, are the same
as in the call to $code optimize_fixed$$ that corresponds to $icode solution$$.

$head information_rcv$$
The return value has prototype
$codei%
	CppAD::mixed::sparse_rcv %information_rcv%
%$$
see $cref/sparse_rcv/typedef/Sparse Types/sparse_rcv/$$.
This is a sparse matrix representation for the
lower triangle of the observed information matrix,
which is symmetric and hence determined by its lower triangle.
Absolute value terms in the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
for the $cref fix_likelihood$$ are not include in this Hessian
because they do not have a derivative (let alone Hessian) at zero.

$children%
	example/user/information_mat.cpp
%$$

$head Example$$
The file $cref information_mat.cpp$$ contains an example and
test of this routine. It returns true for success and false for failure.

$end
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>
# include <cppad/mixed/exception.hpp>

CppAD::mixed::sparse_rcv cppad_mixed::try_information_mat(
	const CppAD::mixed::fixed_solution&  solution             ,
	const d_vector&                      random_opt           )
{	using Eigen::Dynamic;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
	typedef eigen_sparse::InnerIterator                       sparse_itr;
	// solution
	assert( solution.fixed_opt.size()   == n_fixed_ );
	assert( solution.fixed_lag.size()   == n_fixed_ );
	// random_opt
	assert( random_opt.size() == n_random_ );
	//
	// optimal fixed effects
	const d_vector& fixed_opt( solution.fixed_opt );
	//
	// ----------------------------------------------------------------------
	// compute Hessian w.r.t. fixed effect of random part of objective
	eigen_sparse ran_hes;
	if( n_random_ == 0 )
		ran_hes.resize( int(n_fixed_), int(n_fixed_) );
	else
	{	// ------------------------------------------------------------------
		// update the cholesky factor for this fixed and random effect
		update_factor(fixed_opt, random_opt);
		// ------------------------------------------------------------------
		// If Quasi-Newton method was used, must initilaize routines
		// that are only used for the Hessian calculation; see initilaize.cpp
		if( ! init_laplace_obj_done_ )
		{	assert( quasi_fixed_ );
			assert( ! init_laplace_obj_hes_done_ );
			//
			// laplace_obj_fun_
			assert( ! init_laplace_obj_done_ );
			init_laplace_obj(fixed_opt, random_opt);
			assert( init_laplace_obj_done_ );
			//
			// laplace_obj_hes_
			assert( ! init_laplace_obj_hes_done_ );
			init_laplace_obj_hes(fixed_opt, random_opt);
			assert( init_laplace_obj_hes_done_ );
		}
		assert( init_laplace_obj_done_ );
		assert( init_laplace_obj_hes_done_ );
		// ------------------------------------------------------------------
		// Lower triangle of Hessian w.r.t. fixed effects
		// for random part of objective,no constraints
		d_vector w_ran(A_rcv_.nr() + 1);
		w_ran[0] = 1.0;
		for(size_t j = 1; j <=A_rcv_.nr(); j++)
			w_ran[j] = 0.0;
		//
		CppAD::mixed::sparse_mat_info ran_info;
		laplace_obj_hes(
		fixed_opt, random_opt, w_ran, ran_info.row, ran_info.col, ran_info.val
		);
		CppAD::mixed::triple2eigen(
			ran_hes       ,
			n_fixed_      ,
			n_fixed_      ,
			ran_info.row  ,
			ran_info.col  ,
			ran_info.val
		);
	}
	// --------------------------------------------------------------------
	// Lower triangle of Hessian of the fixed likelihood
	eigen_sparse fix_hes;
	size_t n_fix_like = 0;
	if( fix_like_fun_.size_var() == 0 )
		fix_hes.resize( int(n_fixed_) , int(n_fixed_) );
	else
	{	n_fix_like = fix_like_fun_.Range();
		d_vector w_fix(n_fix_like);
		w_fix[0] = 1.0;
		for(size_t j = 1; j < n_fix_like; j++)
			w_fix[0] = 0.0;
		//
		CppAD::mixed::sparse_mat_info fix_info;
		fix_like_hes(
				fixed_opt, w_fix, fix_info.row, fix_info.col, fix_info.val
		);
		CppAD::mixed::triple2eigen(
			fix_hes       ,
			n_fixed_      ,
			n_fixed_      ,
			fix_info.row  ,
			fix_info.col  ,
			fix_info.val
		);
	}
	// -----------------------------------------------------------------------
	// Hessian of total objective (observed information matrix)
	eigen_sparse total_hes = ran_hes + fix_hes;
	//
	// convert from eigen sparse matrix to sparse_rcv
	size_t nr  = total_hes.rows();
	size_t nc  = total_hes.cols();
	size_t nnz = total_hes.nonZeros();
	sparse_rc pattern(nr, nc, nnz);
	d_vector val(nnz);
	size_t k = 0;
	for(size_t j = 0; j < n_fixed_; j++)
	{	for(sparse_itr itr(total_hes, int(j) ); itr; ++itr)
		{	assert( size_t( itr.col() ) == j );
			size_t i = itr.row();
			pattern.set(k, i, j);
			val[k] = itr.value();
			++k;
		}
	}
	assert( k == nnz );
	sparse_rcv total_rcv(pattern);
	for(k = 0;  k < nnz; k++)
		total_rcv.set(k, val[k]);
	//
	return total_rcv;
}
// ---------------------------------------------------------------------------
CppAD::mixed::sparse_rcv cppad_mixed::information_mat(
	const CppAD::mixed::fixed_solution&  solution             ,
	const d_vector&                      random_opt           )
{	sparse_rcv result;
	try
	{	result = try_information_mat(solution, random_opt);
	}
	catch(const CppAD::mixed::exception& e)
	{	std::string error_message = e.message("information_mat");
		fatal_error(error_message);
		assert(false);
	}
	//
	return result;
}
