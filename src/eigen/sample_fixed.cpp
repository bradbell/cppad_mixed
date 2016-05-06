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
$begin sample_fixed$$
$spell
	Cholesky
	covariance
	pairwise
	CppAD
	cppad
	const
	gsl_rng
$$

$section Simulating the Posterior Distribution for Fixed Effects$$

$head Syntax$$
$codei%mixed_object%.sample_fixed(
	%sample%,
	%information_info%,
	%solution%,
	%fixed_lower%,
	%fixed_upper%,
	%random_opt%
)%$$

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head Purpose$$
This routine draw samples from
the asymptotic posterior distribution for the
optimal fixed effects (given the model and the data).

$head manage_gsl_rng$$
It is assumed that
$cref/get_gsl_rng/manage_gsl_rng/get_gsl_rng/$$ will return
a pointer to a GSL random number generator.


$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head sample$$
This argument has prototype
$codei%
	CppAD::vector<double>& %sample%
%$$
and its size is a multiple of
$cref/n_fixed/derived_ctor/n_fixed/$$.
The input value of its elements does not matter.
We define
$codei%
	%n_sample% = %sample_size% / %n_fixed%
%$$
Upon return,
for $codei%i% = 0 , %...%, %n_sample%-1%$$,
$codei%j% = 0 , %...%, %n_fixed%-1%$$,
$codei%
	%sample%[ %i% * %n_fixed% + %j% ]
%$$
is the $th j$$ component of the $th i$$ sample of the
optimal fixed effects $latex \hat{\theta}$$.
These samples are independent for different $latex i$$,
and for fixed $latex i$$, they have the
$cref/constrained covariance/sample_fixed/Theory/Constrained Covariance/$$
$latex D$$.

$head information_info$$
The argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info %information_info%
%$$
This is a sparse matrix representation for the
lower triangle of the observed information matrix corresponding to
$icode solution$$; i.e., the matrix returned by
$codei%
%information_info% = %mixed_object%.information_mat(
	%solution%, %random_options%, %random_lower%, %random_upper%, %random_in%
)%$$

$head solution$$
is the $cref/solution/optimize_fixed/solution/$$
for a the call to $cref optimize_fixed$$ corresponding to
$icode information_info$$.

$head fixed_lower$$
is the same as
$cref/fixed_lower/optimize_fixed/fixed_lower/$$
in the call to $code optimize_fixed$$ that corresponds to $icode solution$$.

$head fixed_upper$$
is the same as
$cref/fixed_upper/optimize_fixed/fixed_upper/$$
in the call to $code optimize_fixed$$ that corresponds to $icode solution$$.

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

$head Theory$$

$subhead Notation$$
Given two random vectors $latex u$$ and $latex v$$,
we use the notation $latex \B{C}( u , v )$$
for the corresponding covariance matrix;
i.e.,
$latex \[
	\B{C}( u , v )
	=
	\B{E} \left( [ u - \B{E} (u) ] [ v - \B{E} (v) ]^\R{T} \right)
\] $$

$subhead Unconstrained Covariance$$
We use $latex \tilde{L} ( \theta )$$ to denote the
$cref/total objective/theory/Objective/Total Objective, L(theta)/$$
with absolute values terms in $cref fix_likelihood$$ excluded.
We use $latex \tilde{\theta}$$ for the unconstrained optimal estimate
of the fixed effects (without absolute value terms) and
approximation its auto-covariance by
$latex \[
	\B{C} ( \tilde{\theta} , \tilde{\theta} )
	=
	H^{-1}
\]$$
where $latex H$$ is the Hessian
corresponding to $icode information_info$$; i.e.,
the observed information matrix.

$subhead Constraints Equations$$
Let $latex n$$ be the number of fixed effects,
$latex m$$ the number of active constraints,
and the equations $latex e( \theta ) = b$$ the active constraints.
Here $latex e : \B{R}^n \rightarrow \B{R}^m$$ and $latex b \in \B{R}^m$$
and the inequality constraints have been converted to equalities at the
active bounds (including the bounds on the fixed effects).
Define the random variable the approximation for $latex e( \theta )$$ by
$latex \[
\tilde{e} ( \theta ) =
e( \hat{\theta} ) + e^{(1)} \left( \hat{\theta} \right)
	\left( \theta - \hat{\theta} \right)
\] $$

$subhead Constrained Covariance$$
We approximate the distribution for
$latex \tilde{\theta}$$ normal,
and the distribution for $latex \hat{\theta}$$
as the conditional distribution of $latex \tilde{\theta}$$ given
the value of $latex \tilde{e} ( \tilde{\theta} )$$; i.e.,
$latex \[
	\B{C} \left( \hat{\theta} \W{,} \hat{\theta} \right)
	=
	\B{C} \left( \tilde{\theta} \W{,} \tilde{\theta} \right)
	-
	\B{C} \left( \tilde{\theta} \W{,} \tilde{e} \right)
	\B{C} \left( \tilde{e}  \W{,} \tilde{e} \right)^{-1}
	\B{C} \left( \tilde{e}  \W{,} \tilde{\theta} \right)
\] $$
Using the notation
$latex D = \B{C} \left( \hat{\theta} \W{,} \hat{\theta} \right)$$,
$latex C = \B{C} \left( \tilde{\theta} \W{,} \tilde{\theta} \right)$$,
$latex E = e^{(1)} \left( \hat{\theta} \right)$$,
we have
$latex \[
	D = C - C E^\R{T} \left( E C E^\R{T} \right)^{-1}  E C
\] $$

$children%example/user/sample_fixed_xam.cpp
%$$
$head Example$$
The file $cref sample_fixed_xam.cpp$$ is an example
and test of $code sample_fixed$$.

$end
------------------------------------------------------------------------------
*/
# include <Eigen/Core>
# include <Eigen/LU>
# include <Eigen/Cholesky>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>

# define DEBUG_PRINT 0

namespace {
	using Eigen::Dynamic;
	using CppAD::mixed::get_gsl_rng;
	typedef Eigen::Matrix<double, Dynamic, Dynamic>     double_mat;
	typedef Eigen::Matrix<double, Dynamic, 1>           double_vec;
	typedef Eigen::Matrix<size_t, Dynamic, 1>           size_vec;
	typedef Eigen::LDLT<double_mat, Eigen::Lower>       double_cholesky;
	typedef Eigen::PermutationMatrix<Dynamic, Dynamic, int>  permutation_mat;
	//
# if DEBUG_PRINT
	void print(const char* name , const double_mat& mat)
	{	std::cout << "\n" << name << " =\n" << mat << "\n"; }
	void print(const char* name , double_vec& vec)
	{	std::cout << "\n" << name << "^T = " << vec.transpose() << "\n"; }
	void print(const char* name , size_vec& vec)
	{	std::cout << "\n" << name << "^T = " << vec.transpose() << "\n"; }
# endif
}


void cppad_mixed::sample_fixed(
	d_vector&                            sample               ,
	const CppAD::mixed::sparse_mat_info& information_info     ,
	const CppAD::mixed::fixed_solution&  solution             ,
	const d_vector&                      fixed_lower          ,
	const d_vector&                      fixed_upper          ,
	const d_vector&                      random_opt           )
{
	// number of fixed constraints
	size_t n_fix_con = 0;
	if( fix_con_fun_.size_var() != 0 )
		n_fix_con = fix_con_fun_.Range();
	//
	// sample
	assert( sample.size() > 0 );
	assert( sample.size() % n_fixed_ == 0 );
	// information_info
	assert( information_info.row.size() == information_info.col.size() );
	assert( information_info.row.size() == information_info.val.size() );
	// solution
	assert( solution.fixed_opt.size() == n_fixed_ );
	assert( solution.fixed_lag.size() == n_fixed_ );
	assert( solution.fix_con_lag.size() == n_fix_con );
	assert( solution.ran_con_lag.size() == n_ran_con_ );
	// random_opt
	assert( random_opt.size() == n_random_ );
	//
	// number of samples
	size_t n_sample = sample.size() / n_fixed_;
	//
	// optimal fixed effects
	const d_vector& fixed_opt( solution.fixed_opt );
	// -----------------------------------------------------------------------
	// update the cholesky factor for this fixed and random effect
	if( n_random_ > 0 )
		update_factor(fixed_opt, random_opt);
	// -----------------------------------------------------------------------
	// Create con_mat
	//
	// number of active bound constraints
	size_t n_bnd_active = 0;
	for(size_t j = 0; j < n_fixed_; j++)
	{	if( solution.fixed_lag[j] != 0.0 )
			n_bnd_active++;
	}
	// number fixed constraints active
	size_t n_fix_active = 0;
	size_vec fix_active_index(n_fix_con);
	for(size_t i = 0; i < n_fix_con; i++)
	{	fix_active_index[i] = n_fixed_;
		if( solution.fix_con_lag[i] != 0.0 )
			fix_active_index[i] = n_fix_active++;
	}
	// number of random constraints active
	size_t n_ran_active = n_ran_con_;
	//
	// matrix with all the active constraints
	size_t n_con_active = n_bnd_active  + n_fix_active + n_ran_active;
	double_mat con_mat = double_mat::Zero(n_con_active, n_fixed_);
	size_t con_row = 0;
	//
	// put the bounds in con_mat
	for(size_t j = 0; j < n_fixed_; j++) if( solution.fixed_lag[j] != 0.0 )
	{	con_mat(con_row, j) = 1.0;
		++con_row;
	}
	// put fixed constraints in con_mat
	if( n_fix_con > 0 )
	{	assert( con_row == n_bnd_active );
		//
		// jacobian of the fixed constraints
		CppAD::mixed::sparse_mat_info fix_con_info;
		fix_con_jac(
			fixed_opt, fix_con_info.row, fix_con_info.col, fix_con_info.val
		);
		size_t K = fix_con_info.row.size();
		for(size_t k = 0; k < K; k++)
		{	size_t r = fix_con_info.row[k];
			if( solution.fix_con_lag[r] != 0.0 )
			{	assert( fix_active_index[r] != n_fixed_ );
				size_t i = fix_active_index[r];
				size_t c  = fix_con_info.col[k];
				double v  = fix_con_info.val[k];
				con_mat(con_row + i, c) = v;
			}
		}
		con_row += n_fix_active;
	}
	// put random constraints in con_mat
	CppAD::mixed::sparse_mat_info ran_con_info;
	if( n_ran_active > 0 )
	{	assert( con_row == n_bnd_active + n_fix_active );
		//
		// jacobian of the random constraints
		// sparsity pattern
		ran_con_jac(fixed_opt, random_opt, ran_con_info);
		// values
		ran_con_jac(fixed_opt, random_opt, ran_con_info);
		//
		size_t K = ran_con_info.row.size();
		for(size_t k = 0; k < K; k++)
		{	size_t r = ran_con_info.row[k];
			size_t c  = ran_con_info.col[k];
			double v  = ran_con_info.val[k];
			con_mat(con_row + r, c) = v;
		}
	}
	// -----------------------------------------------------------------------
	// compute conditional covariance
	// ----------------------------------------------------------------------
	//
	// information matrix
	double_mat info_mat = double_mat::Zero(n_fixed_, n_fixed_);
	for(size_t k = 0; k < information_info.row.size(); k++)
	{	// note only lower triangle is stored in information_info
		size_t r = information_info.row[k];
		size_t c = information_info.col[k];
		double v = information_info.val[k];
		info_mat(r, c) = v;
		info_mat(c, r) = v;
	}
	//
	// covariance matrix with out constraints
	double_mat C = double_mat( info_mat.inverse() );
	//
	// conditional covariance
	double_mat& E(con_mat);
	double_mat  EC    = E * C;
	double_mat  ECET  = EC * E.transpose();
	double_mat  D     = C - EC.transpose() * ECET.inverse() * EC;
	//
	// LDLT factorizaton of D
	double_cholesky cholesky;
	cholesky.compute(D);
	//
	// diagonal elements of LDLT factorization
	double_vec diag      = cholesky.vectorD();
	double_vec diag_root(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
	{	if( diag[j] > 0.0 )
			diag_root[j] = std::sqrt( diag[j] );
		else
			diag_root[j] = 0.0;
	}
	double_mat L      = cholesky.matrixL();
	permutation_mat P = permutation_mat( cholesky.transpositionsP() );
	// -----------------------------------------------------------------------
	// Simulate the samples
	// -----------------------------------------------------------------------
	for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
	{	double_vec w(n_fixed_);
		// simulate a normal with mean zero and variance sqrt{D(j,j)}
		for(size_t j = 0; j < n_fixed_; j++)
			w[j] = diag_root[j] * gsl_ran_gaussian(get_gsl_rng(), 1.0);
		// multily by Cholesky factor
		double_vec s = P.transpose() * L * w;
		//
		for(size_t j = 0; j < n_fixed_; j++)
		{
			// store corresponding sample
			sample[ i_sample * n_fixed_ + j] = fixed_opt[j] + s[j];
			//
			// a constraint that is not active might go out of bounds
			s[j] = std::min(s[j], fixed_upper[j]);
			s[j] = std::max(s[j], fixed_lower[j]);
		}
	}
	// -----------------------------------------------------------------------
	return;
}
