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
$$

$section Simulating the Posterior Distribution for the Fixed Effects$$

$head Syntax$$
$icode%correlation% = %mixed_object%.sample_fixed(
	%sample%,
	%non_zero%,
	%information_info%,
	%solution%,
	%fixed_lower%,
	%fixed_upper%,
	%fixed_constraint_lower%,
	%fixed_constraint_upper%,
	%random_lower%,
	%random_upper%,
	%random_in%,
)%$$

$head Purpose$$
This is just an example of how one might draw samples from
the asymptotic posterior distribution for the
optimal fixed effects (given the model and the data).
It is not implemented.

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

$head sample$$
This argument has prototype
$codei%
	CppAD::vector<double>& %sample%
%$$
and its size is a multiple of
$cref/n_fixed/derived_ctor/n_fixed/$$.
We define
$codei%
	%n_sample% = %sample_size% / %n_fixed%
%$$
The input value of its elements does not matter.
Upon return,
for $codei%i% = 0 , %...%, %n_sample%-1%$$,
$codei%j% = 0 , %...%, %n_fixed%-1%$$,
$codei%
	%sample%[ %i% * %n_fixed% + %j% ]
%$$
is the $th j$$ component of the $th i$$ sample of the
optimal fixed effects $latex \hat{\theta}$$.

$head non_zero$$
This argument has prototype
$codei%
	double  %non_zero%
%$$
and is a value between zero and one.
It specifies the fraction of the elements in
the posterior covariance that are included in the simulation.
Diagonal elements are always included in the simulation, so
$codei%
	%non_zero% >= 1.0 - %n_fixed% / (%n_fixed% * %n_fixed%)
%$$
has the same effects at $icode non_zero$$ equal to one.

$head information_info$$
The argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info %information_info%
%$$
This is a sparse matrix representation for the
lower triangle of the observed information matrix; see
$cref information_mat$$.

$head solution$$
is the $cref/solution/optimize_fixed/solution/$$
for a previous call to $cref optimize_fixed$$.

$head random_options$$
is the $cref/random_options/optimize_fixed/random_options/$$
for a previous call to $code optimize_fixed$$.

$head fixed_lower$$
is the $cref/fixed_lower/optimize_fixed/fixed_lower/$$
for a previous call to $code optimize_fixed$$.

$head fixed_upper$$
is the $cref/fixed_upper/optimize_fixed/fixed_upper/$$
for a previous call to $code optimize_fixed$$.

$head fix_constraint_lower$$
is the $cref/fix_constraint_lower/optimize_fixed/fix_constraint_lower/$$
for a previous call to $code optimize_fixed$$.

$head fix_constraint_upper$$
is the $cref/fix_constraint_upper/optimize_fixed/fix_constraint_upper/$$
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

$head correlation$$
The return value has prototype
$codei%
	double %correlation%
%$$
This is the largest absolute correlation
that is converted to zero before simulating the posterior samples.
If $icode%non_zero% = 1.0%$$, no values are treated as zero
and $icode%correlation% = 0.0%$$.


$head Unconstrained Covariance$$
We use $latex L ( \theta )$$ to denote the
$cref/total objective/theory/Objective/Total Objective, L(theta)/$$.
If $latex \hat{\theta}$$ is the optimal value for the fixed effects,
the corresponding estimate for the covariance of $latex \hat{\theta}$$ is
given by
$latex \[
	\B{C} ( \hat{\theta}, \hat{\theta} )
	=
	L^{(2)} ( \hat{\theta} )^{-1}
\]$$
Absolute value terms in the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
for the $cref fix_likelihood$$ are not include in this Hessian
(because they do not have a derivative, let alone Hessian, at zero).


$head Approximate Constraints$$
Let $latex n$$ be the number of fixed constraints,
$latex n$$ the number of active constraints,
and the equations $latex e( \theta ) = b$$ the active constraints
where $latex e : \B{R}^n \rightarrow \B{R}^m$$, $latex b \in \B{R}^m$$,
and the inequality constraints have been converted to equalities at the
active bounds.
Define the random variable $latex \tilde{e}$$ as the affine
approximation for $latex e( \theta )$$.
$latex \[
	\tilde{e} = b + e^{(1)} ( \hat{\theta} ) ( \theta - \hat{\theta} )
\] $$

$head Constrained Covariance$$
We approximate the distribution for $latex \hat{\theta}$$ as normal,
hence $latex \tilde{e}$$ is also normal.
Furthermore, we approximate the distribution for $latex \hat{\theta}$$,
by is conditional distribution given that $latex \tilde{e}$$
is equal to $latex b$$.
In other words,
$latex \[
	\B{C} ( \hat{\theta} | b \W{,} \hat{\theta} | b )
	=
	\B{C} ( \hat{\theta} \W{,} \hat{\theta} )
	-
	\B{C} ( \hat{\theta} \W{,} \tilde{e} )
	\B{C} ( \tilde{e}  \W{,} \tilde{e} )^{-1}
	\B{C} ( \tilde{e}  \W{,} \hat{\theta} )
\] $$
Using the notation
$latex C = L^{(2)} ( \hat{\theta} )^{-1}$$
and
$latex E = e^{(1)} ( \hat{\theta} )$$,
we define
$latex \[
	D
	=
	\B{C} ( \hat{\theta} | b \W{,} \hat{\theta} | b )
	=
	C - C E^\R{T} \left( E C E^\R{T} \right)^{-1}  E C
\] $$
We use $latex D( \alpha )$$ for the matrix where the
off diagonal elements of $latex D$$ corresponding to an absolute correlation
less than $latex \alpha$$ are replaced by zero.
To be specific,
$latex \[
D_{i,j} ( \alpha ) = \left\{ \begin{array}{ll}
0 & \R{if} \; i \neq j \; \R{and} \; \alpha \geq | D_{i,j} | / \sqrt{D_ii D_jj}
\\
D_{i,j}       & \R{otherwise}
\end{array} \right.
\] $$
We define $latex \bar{D} = D ( \alpha )$$ where $latex \alpha$$
is smallest value, greater than or equal zero,
such that the fraction of non-zero values in
$latex \bar{D}$$ is greater than or equal $icode non_zero$$.
The return value $cref/correlation/sample_fixed/correlation/$$ is
the maximum value of $latex | D_{i,j} | / \sqrt{ D_ii D_jj }$$
that is converted to zero in the definition of $latex \bar{D}$$.

$head Simulation$$
We use $latex \sqrt{\bar{D}}$$ to denote a Cholesky factor of
$latex \bar{D}$$; i.e.,
$latex \bar{D} = \sqrt{\bar{D}} \sqrt{\bar{D}}^\R{T}$$.
Further suppose that for $icode%i% = 0 , %...%, %n_sample%-1%$$,
$latex w^i \in \B{R}^n$$ is simulated as independent normal random vectors
with identity covariance matrix.
The $th i$$ the samples for the fixed effects estimate $latex \hat{\theta}$$
is $latex \sqrt{\bar{D}} w^i$$.

$end
------------------------------------------------------------------------------
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/triple2eigen.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>

double cppad_mixed::sample_fixed(
	d_vector&                            sample               ,
	double                               non_zero             ,
	const CppAD::mxied::sparse_mat_info& information_info     ,
	const CppAD::mixed::fixed_solution&  solution             ,
	const std::string&                   random_options       ,
	const d_vector&                      fixed_lower          ,
	const d_vector&                      fixed_upper          ,
	const d_vector&                      fix_constraint_lower ,
	const d_vector&                      fix_constraint_upper ,
	const d_vector&                      random_lower         ,
	const d_vector&                      random_upper         ,
	const d_vector&                      random_in            )
{	using Eigen::Dynamic;
	using CppAD::mixed::get_gsl_rng;
	typedef Eigen::Matrix<double, Dynamic, Dynamic>           eigen_matrix;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
	typedef Eigen::SimplicialLLT<eigen_sparse, Eigen::Lower>  eigen_cholesky;
	typedef eigen_sparse::InnerIterator                       sparse_itr;
	//
	// number of fixed constraints
	size_t n_fix_con = 0;
	if( fix_con_fun_.size_var() != 0 )
		n_fix_con = fix_con_fun_.Range();
	//
	// sample
	assert( sample.size() > 0 );
	assert( sample.size() % n_fixed_ == 0 );
	// non_zero
	assert( 0.0 <= non_zero && non_zero <= 1.0 );
	// solution
	assert( solution.fixed_opt.size() == n_fixed_ );
	assert( solution.fixed_lag.size() == n_fixed_ );
	assert( solution.fix_con_lag.size() == n_fix_con );
	assert( solution.ran_con_lag.size() == n_ran_con_ );
	// fixed_(lower and upper)
	assert( fixed_lower.size() == n_fixed_ );
	assert( fixed_upper.size() == n_fixed_ );
	// fix_constraint(lower and upper)
	assert( fix_constraint_lower.size() == n_fix_con );
	assert( fix_constraint_upper.size() == n_fix_con );
	// random_(lower, upper, in)
	assert( random_lower.size() == n_random_ );
	assert( random_upper.size() == n_random_ );
	assert( random_in.size() == n_random_ );
	//
	// number of samples
	size_t n_sample = sample.size() / n_fixed_;
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
		init_newton_atom_done_ = true;
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
	// full covariance (inverse of the observed information matrix)
	eigen_sparse total_hes = CppAD::mixed::triple2eigen(
		n_fixed_              ,
		n_fixed_              ,
		information_info.row  ,
		information_info.col  ,
		information_info.val
	);
	//
	// identity matrix
	eigen_sparse eye(n_fixed_, n_fixed_);
	for(size_t i = 0; i < n_fixed_; i++)
		eye.insert(i, i) = 1.0;
	//
	// Inverse of total_hes is our approximate unconstrained covariance
	eigen_cholesky cholesky;
	cholesky.compute(total_hes);
	eigen_sparse full_cov = cholesky.solve(eye);
	// -----------------------------------------------------------------------
	// Subtract fixed and random constraints from the full covariance
	//
	// jacobian of the fixed constraints
	CppAD::mixed::sparse_mat_info fix_con_info;
	if( n_fix_con > 0 )
	{	fix_con_jac(
			fixed_opt, fix_con_info.row, fix_con_info.col, fix_con_info.val
		);
	}
	// jacobian of the random constraints
	CppAD::mixed::sparse_mat_info ran_con_info;
	if( n_ran_con_ > 0 )
	{	ran_con_jac(fixed_opt, random_opt, ran_con_info);
	}
	// sparsity triple for both fixed and random constraints
	CppAD::mixed::sparse_mat_info con_info;
	for(size_t k = 0; k < fix_con_info.row.size(); k++)
	{	con_info.row.push_back( fix_con_info.row[k] );
		con_info.col.push_back( fix_con_info.col[k] );
		con_info.val.push_back( fix_con_info.val[k] );
	}
	for(size_t k = 0; k < ran_con_info.row.size(); k++)
	{	con_info.row.push_back( ran_con_info.row[k] );
		con_info.col.push_back( ran_con_info.col[k] );
		con_info.val.push_back( ran_con_info.val[k] );
	}
	//
	if( con_info.row.size() > 0 )
	{	eigen_sparse E = CppAD::mixed::triple2eigen(
			n_fixed_      ,
			n_fixed_      ,
			con_info.row  ,
			con_info.col  ,
			con_info.val
		);
		eigen_sparse EC = E * full_cov;
		cholesky.compute( EC * E.transpose() );
		eye.resize(n_fix_con, n_fix_con);
		for(size_t i = 0; i < n_fix_con; i++)
			eye.insert(i, i) = 1.0;
		eigen_sparse ECE_inv = cholesky.solve(eye);
		full_cov -= EC.transpose() * ECE_inv * EC;
	}
	// -----------------------------------------------------------------------
	// full2reduced, reduced2full, n_reduced
	CppAD::vector<size_t> full2reduced(n_fixed_), reduced2full(0);
	for(size_t i = 0; i < n_fixed_; i++)
	{	if( solution.fixed_lag[i] != 0.0 )
		{	// this variable is not in the reduced matrix
			full2reduced[i] = n_fixed_;
		}
		else
		{	// mapping from full index to reduced index
			full2reduced[i] = reduced2full.size();
			// mapping from reduced index to full index
			reduced2full.push_back(i);
		}
	}
	size_t n_reduced = reduced2full.size();
	// ----------------------------------------------------------------------
	// determine which components of the reduced covaraince to zero out
	//
	// only really need lower triangle, but use full matrix to simplify code
	size_t n_fixed_sq = n_fixed_ * n_fixed_;
	CppAD::vector<double> diagonal(n_fixed_), ratio(0), row(0), col(0);
	for(size_t j = 0; j < n_fixed_; j++)
	{	diagonal[j]               = 0.0;
		for(sparse_itr itr(full_cov, j); itr; ++itr)
		{	assert( size_t( itr.col() ) == j );
			if( itr.row() == itr.col() )
				diagonal[j] = itr.value();
		}
		assert( diagonal[j] > 0.0 );
	}
	double eps = 100. * std::numeric_limits<double>::epsilon();
	for(size_t j = 0; j < n_fixed_; j++)
	{	for(sparse_itr itr(full_cov, j); itr; ++itr)
		{	assert( size_t( itr.col() ) == j );
			size_t r = itr.row();
			size_t c = itr.col();
			if( r == c )
			{	// diagonal eleemnts all have ratio 1.0
				ratio.push_back( 1.0 );
				row.push_back(r);
				col.push_back(c);
			}
			// off diagonal elements are zero for variables at bounds
			else if( full2reduced[j] != n_fixed_ )
			{	double scale = std::sqrt( diagonal[r] * diagonal[c] );
				ratio.push_back( std::min( itr.value() / scale, 1.0 - eps ) );
				row.push_back(r);
				col.push_back(c);
			}
		}
	}
	CppAD::vector<size_t> ind(n_fixed_sq);
	CppAD::index_sort(ratio, ind);
	//
	double correlation = 0.0;
	size_t n_zero      = (1.0 - non_zero) * n_fixed_sq;
	n_zero             = std::min(n_zero, n_fixed_sq - n_fixed_);
	if( n_fixed_sq - ratio.size() < n_zero )
	{	size_t n_change = n_zero - (n_fixed_sq - ratio.size());
		for(size_t k = 0; k < n_change; k++)
			full_cov.coeffRef( row[ind[k]], col[ind[k]] ) = 0.0;
		// maximum correlation that was changed to zero
		correlation = ratio[ ind[ n_change - 1 ] ];
	}
	// -----------------------------------------------------------------------
	// Bound constrained variables get removed form the covariance
	//
	//
	// reduced_cov
	eigen_sparse reduced_cov(n_reduced, n_reduced);
	for(size_t r_col = 0; r_col < n_reduced; r_col++)
	{	size_t f_col = reduced2full[r_col];
		for(sparse_itr itr(full_cov, f_col); itr; ++itr)
		{	size_t f_row = itr.row();
			size_t r_row = full2reduced[f_row];
			if( r_row != n_fixed_ )
			{	assert( r_row < n_reduced );
				reduced_cov.insert(r_row, r_col);
			}
		}
	}
	// -----------------------------------------------------------------------
	// Simulate the samples
	//
	// Cholesky factor for reduced covariance
	cholesky.compute(reduced_cov);
	//
	for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
	{	eigen_matrix w(n_reduced, 1);
		// simulate a normal with mean zero and variance on
		for(size_t j = 0; j < n_reduced; j++)
			w(j, 0) = gsl_ran_gaussian(get_gsl_rng(), 1.0);
		// multily by Cholesky factor
		eigen_matrix s = cholesky.matrixL() * w;
		// store in sample
		for(size_t j = 0; j < n_reduced; j++)
			sample[ i_sample * n_fixed_ + reduced2full[j] ] = s(j, 0);
		for(size_t j = 0; j < n_fixed_; j++)
		{	if( full2reduced[j] == n_fixed_ )
			{	if( solution.fixed_lag[j] > 0.0 )
					sample[ i_sample * n_fixed_ + j ] = fixed_lower[j];
				else
				{	assert( solution.fixed_lag[j] < 0.0 );
					sample[ i_sample * n_fixed_ + j ] = fixed_upper[j];
				}
			}
		}
	}
	// -----------------------------------------------------------------------
	return correlation;
}
