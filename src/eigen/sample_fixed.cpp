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
$$

$section Simulating the Posterior Distribution for the Fixed Effects$$

$head Syntax$$
$icode%correlation% = %mixed_object%.sample_fixed(
	%sample%,
	%non_zero%,
	%solution%,
	%fixed_lower%,
	%fixed_upper%,
	%fixed_constraint_lower%,
	%fixed_constraint_upper%,
	%random_lower%,
	%random_upper%,
	%random_in%,
)%$$

$head Under Construction$$

$head Purpose$$
Sample the asymptotic posterior distribution for the
optimal fixed effects (given the model and the data).

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
It specifies the fraction of off diagonal elements in
the posterior covariance that are included in the simulation.


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
that is treated as zero when simulating the posterior samples.
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
such that the fraction of non-zero off diagonal values in
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

double cppad_mixed::sample_fixed(
	d_vector&                            sample               ,
	double                               non_zero             ,
	const CppAD::mixed::fixed_solution&  solution             ,
	const std::string&                   random_options       ,
	const d_vector&                      fixed_lower          ,
	const d_vector&                      fixed_upper          ,
	const d_vector&                      fix_constraint_lower ,
	const d_vector&                      fix_constraint_upper ,
	const d_vector&                      random_lower         ,
	const d_vector&                      random_upper         ,
	const d_vector&                      random_in            )
{
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor>      eigen_sparse;
	typedef Eigen::SimplicialLDLT<eigen_sparse, Eigen::Lower> eigen_cholesky;
	typedef eigen_sparse::InnerIterator                       sparse_itr;
	//
	// sample
	assert( sample.size() > 0 );
	assert( sample.size() % n_fixed_ == 0 );
	// non_zero
	assert( 0.0 <= non_zero && non_zero <= 1.0 );
	// solution
	assert( solution.fixed_opt.size() == n_fixed_ );
	assert( solution.fixed_lag.size() == n_fixed_ );
	assert( solution.fix_con_lag.size() == fix_con_fun_.Range() );
	assert( solution.ran_con_lag.size() == n_ran_con_ );
	// fixed_(lower and upper)
	assert( fixed_lower.size() == n_fixed_ );
	assert( fixed_upper.size() == n_fixed_ );
	// fix_constraint(lower and upper)
	assert( fix_constraint_lower.size() == fix_con_fun_.Range() );
	assert( fix_constraint_upper.size() == fix_con_fun_.Range() );
	// random_(lower, upper, in)
	assert( random_lower.size() == n_random_ );
	assert( random_upper.size() == n_random_ );
	assert( random_in.size() == n_random_ );
	//
	// number of samples
	// size_t n_sample = sample.size() / n_fixed_;
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
	// Compute our approximation for unconstrained covariance full_cov
	//
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
	//
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
	fix_con_hes(
			fixed_opt, w_fix, fix_info.row, fix_info.col, fix_info.val
	);
	eigen_sparse fix_hes = CppAD::mixed::triple2eigen(
		n_fixed_      ,
		n_fixed_      ,
		fix_info.row  ,
		fix_info.col  ,
		fix_info.val
	);
	//
	// Hessian of total objective (observed information matrix)
	eigen_sparse total_hes = ran_hes + fix_hes;
	//
	// identity matrix
	eigen_sparse eye(n_fixed_, n_fixed_);
	for(size_t i = 0; i < n_fixed_; i++)
		eye.insert(i, i);
	//
	// Inverse of total_hes is our approximate unconstrained covariance
	eigen_cholesky cholesky;
	cholesky.compute(total_hes);
	eigen_sparse full_cov = cholesky.solve(eye);
	// -----------------------------------------------------------------------
	// Bound constrained variables get removed form the covariance
	//
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
	// under construction
	double sum = 0.0;
	for(sparse_itr itr(reduced_cov, 0); itr; ++itr)
		sum += itr.value();
	return sum;
}
