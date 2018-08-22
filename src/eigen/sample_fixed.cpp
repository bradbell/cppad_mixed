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
$begin sample_fixed$$
$spell
	Cholesky
	covariance
	pairwise
	CppAD
	cppad
	const
	gsl_rng
	rcv
$$

$section Sample Posterior for Fixed Effects$$

$head Syntax$$
$icode%mixed_object%.sample_fixed(
	%sample%,
	%information_rcv%,
	%solution%,
	%fixed_lower%,
	%fixed_upper%,
	%random_opt%
)%$$

$head See Also$$
$cref sample_random$$, $cref sample_conditional$$

$head Prototype$$
$srcfile%src/eigen/sample_fixed.cpp
	%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head Purpose$$
This routine draws samples from
the asymptotic posterior distribution for the
fixed effects (given the model and the data).

$head manage_gsl_rng$$
It is assumed that
$cref/get_gsl_rng/manage_gsl_rng/get_gsl_rng/$$ will return
a pointer to a GSL random number generator.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head sample$$
The size $icode%sample%.size()%$$ is a multiple of
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
$cref/implicit covariance/sample_fixed/Theory/Implicit Covariance/$$.

$head information_rcv$$
This is a sparse matrix representation for the
lower triangle of the observed information matrix corresponding to
$icode solution$$; i.e., the matrix returned by
$codei%
%information_rcv% = %mixed_object%.information_mat(
	%solution%, %random_options%, %random_lower%, %random_upper%, %random_in%
)%$$

$head solution$$
is the $cref/solution/optimize_fixed/solution/$$
for a the call to $cref optimize_fixed$$ corresponding to
$icode information_rcv$$.

$head fixed_lower$$
is the same as
$cref/fixed_lower/optimize_fixed/fixed_lower/$$
in the call to $code optimize_fixed$$ that corresponding to $icode solution$$.

$head fixed_upper$$
is the same as
$cref/fixed_upper/optimize_fixed/fixed_upper/$$
in the call to $code optimize_fixed$$ that corresponding to $icode solution$$.

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

$head Positive Definite$$
If the $cref/implicit information/sample_fixed/Theory/Implicit Information/$$
matrix is not positive definite, at least one fixed effects is not determined
by the model plus the data.
You should be able to recognize it because its
upper and lower limits are not equal, but the limits appear often
in the simulated values (or it is uniformly distributed between these limits).

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

$subhead Unconstrained Subset Covariance$$
We use $latex \tilde{L} ( \theta )$$ to denote the
$cref/total objective/theory/Objective/Total Objective, L(theta)/$$
as a function of $latex \theta$$ and
where the absolute values terms in $cref fix_likelihood$$ are excluded.
We use $latex \tilde{\theta}$$ for the unconstrained optimal estimate
of the subset of fixed effects and
approximate its auto-covariance by
$latex \[
	\B{C} ( \tilde{\theta} , \tilde{\theta} )
	=
	H^{-1}
\]$$
Here $latex H$$ is the Hessian corresponding to
the observed information matrix $icode information_rcv$$.

$subhead Approximate Constraint Equations$$
Let $latex n$$ be the number of fixed effects in $latex \theta$$,
$latex m$$ the number of active constraints (not counting bounds),
and the equations $latex e( \theta ) = b$$ be the active constraints.
Here $latex e : \B{R}^n \rightarrow \B{R}^m$$ and $latex b \in \B{R}^m$$
and the inequality constraints have been converted to equalities at the
active bounds.
Indices $icode j$$ for which $icode%fixed_lower%[%j%]%$$
is equal to $icode%fixed_upper%[%j%]%$$ are considered even of the
corresponding Lagrange multiplier is zero.
Define the approximate constraint equation by
$latex \[
b =
e \left( \hat{\theta} \right) + e^{(1)} \left( \hat{\theta} \right)
	\left( \theta - \hat{\theta} \right)
\] $$
where $latex \hat{\theta}$$ is the optional estimate
for the fixed effects $icode%solution%.fixed_opt%$$.

$subhead Implicit Information$$
We apply the implicit function theorem to the approximate constraint
equation above to get a representation:
$latex \[
	\theta_D = C \theta_I + d
\] $$
Here $latex D \in \B{Z}_+^m$$ is the dependent subset of $latex \theta$$,
$latex I \in \B{Z}_+^{n-m}$$ is the independent subset of $latex \theta$$,
$latex D \cap I = \emptyset$$,
$latex C \in \B{R}^{m \times (n - m)}$$,
$latex d \in \B{R}^m$$.
For any $latex \theta$$ satisfying the approximate constraint equation,
the corresponding components $latex \theta_D$$ and $latex \theta_I$$
satisfy the equation above.
The maximum likelihood problem corresponding to the approximate constraints,
as a function of $latex \theta_I$$ alone,
would have the following observed information matrix:
$latex \[
	H_{I,I} + H_{I,D} C + C^\R{T} H_{D,I} + C^\R{T} H_{D,D} C
\] $$
This matrix is the observed implicit information matrix.

$subhead Implicit Covariance$$
The inverse of the implicit information matrix
is called the implicit covariance.

$children%example/user/sample_fixed.cpp
	%src/eigen/sample_conditional.cpp
%$$
$head Example$$
The file $cref sample_fixed.cpp$$ is an example
and test of $code sample_fixed$$.

$head Other Method$$
The routine $cref sample_conditional$$, is a old method that is no
longer used for computing these samples.

$head 2DO$$
This routine uses dense matrices, perhaps it would be useful
to converting this (and $cref undetermined$$) to all sparse matrices.

$end
------------------------------------------------------------------------------
*/
# include <Eigen/Core>
# include <Eigen/LU>
# include <Eigen/Cholesky>
# include <gsl/gsl_randist.h>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <cppad/mixed/undetermined.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>

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
// -------------------------------------------------------------------------
void cppad_mixed::try_sample_fixed(
	CppAD::vector<double>&                 sample               ,
	const d_sparse_rcv&                    information_rcv      ,
	const CppAD::mixed::fixed_solution&    solution             ,
	const CppAD::vector<double>&           fixed_lower          ,
	const CppAD::vector<double>&           fixed_upper          ,
	const CppAD::vector<double>&           random_opt           )
{
	// number of fixed constraints
	size_t n_fix_con = 0;
	if( fix_con_fun_.size_var() != 0 )
		n_fix_con = fix_con_fun_.Range();
	//
	// sample
	assert( sample.size() > 0 );
	assert( sample.size() % n_fixed_ == 0 );
	// solution
	assert( solution.fixed_opt.size() == n_fixed_ );
	assert( solution.fixed_lag.size() == n_fixed_ );
	assert( solution.fix_con_lag.size() == n_fix_con );
	assert( solution.ran_con_lag.size() == A_rcv_.nr() );
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
	// Determine the subset of variables that do not have active bounds
	// mapping from fixed index to subset index and back
	CppAD::vector<size_t> fixed2subset(n_fixed_);
	size_t n_subset = 0;
	for(size_t j = 0; j < n_fixed_; j++)
	{	if( solution.fixed_lag[j] != 0.0 )
			fixed2subset[j] = n_fixed_;
		else if( fixed_lower[j] == fixed_upper[j] )
			fixed2subset[j] = n_fixed_;
		else
			fixed2subset[j] = n_subset++;
	}
	assert( n_subset <= n_fixed_ );
	// -----------------------------------------------------------------------
	// Represent approximate constraint as
	// e = E * alpha
	// where alpha is subset of fixed effecst that do not have active bounds
	//
	// number fixed constraints active
	size_t n_fix_active = 0;
	size_vec fix_active_index(n_fix_con);
	for(size_t i = 0; i < n_fix_con; i++)
	{	fix_active_index[i] = n_fixed_;
		if( solution.fix_con_lag[i] != 0.0 )
			fix_active_index[i] = n_fix_active++;
	}
	// random constraints are always active
	size_t n_ran_active = A_rcv_.nr();
	//
	// This A is the matrix with all the active constraints.
	// This includes fixed and random, but not bound, constraints.
	size_t n_con_active = n_fix_active + n_ran_active;
	double_mat A = double_mat::Zero(n_con_active, n_subset);
	double_vec b = double_vec::Zero(n_con_active);
	size_t con_row = 0;
	//
	// put fixed constraints in A
	if( n_fix_con > 0 )
	{	// righth and side for fixed constraints
		d_vector vec  = fix_con_eval(fixed_opt);
		size_t i_active = 0;
		for(size_t r = 0; r < n_fix_con; r++)
		{	if( solution.fix_con_lag[r] != 0.0 )
				b[i_active++] = vec[r];
		}
		// jacobian of the fixed constraints
		CppAD::mixed::sparse_mat_info fix_con_info;
		fix_con_jac(
			fixed_opt, fix_con_info.row, fix_con_info.col, fix_con_info.val
		);
		size_t K = fix_con_info.row.size();
		for(size_t k = 0; k < K; k++)
		{	size_t r  = fix_con_info.row[k];
			size_t c  = fix_con_info.col[k];
			double v  = fix_con_info.val[k];
			//
			bool r_active  = solution.fix_con_lag[r] != 0.0;
			bool in_subset = fixed2subset[c] != n_fixed_;
			if( r_active & in_subset )
			{	assert( fix_active_index[r] != n_fixed_ );
				size_t i  = fix_active_index[r];
				size_t j  = fixed2subset[c];
				A(con_row + i, j) = v;
			}
		}
		con_row += n_fix_active;
	}
	// put random constraints in A
	// (right hand side b for random constraints is zero)
	d_sparse_rcv ran_con_rcv;
	if( n_ran_active > 0 )
	{	assert( con_row == n_fix_active );
		//
		// jacobian of the random constraints
		// sparsity pattern
		ran_con_jac(fixed_opt, random_opt, ran_con_rcv);
		// values
		ran_con_jac(fixed_opt, random_opt, ran_con_rcv);
		//
		size_t K = ran_con_rcv.nnz();
		for(size_t k = 0; k < K; k++)
		{	size_t r       = ran_con_rcv.row()[k];
			size_t c       = ran_con_rcv.col()[k];
			double v       = ran_con_rcv.val()[k];
			bool in_subset = fixed2subset[c] != n_fixed_;
			if( in_subset )
			{	size_t j  = fixed2subset[c];
				A(con_row + r, j) = v;
			}
		}
	}
	double tol = 1e3 * std::numeric_limits<double>::epsilon();
	size_t nD = n_con_active;
	size_t nI = n_subset - n_con_active;
	size_vec   D(nD), I(nI);
	double_mat C(nD, nI);
	double_vec e(nD);
	CppAD::mixed::undetermined(A, b, tol, D, I, C, e);
	// -----------------------------------------------------------------------
	// map from subset index to dependent and independent variable indices
	CppAD::vector<size_t> subset2independent(n_subset);
	CppAD::vector<size_t> subset2dependent(n_subset);
	for(size_t i = 0; i < nI; i++)
	{	subset2independent[ I[i] ] = i;
		subset2dependent[ I[i] ]   = nD;
	}
	for(size_t i = 0; i < nD; i++)
	{	subset2independent[ D[i] ] = nI;
		subset2dependent[ D[i] ]   = i;
	}
	// -----------------------------------------------------------------------
	// compute implicit information
	double_mat H_II = double_mat::Zero(nI, nI);
	double_mat H_ID = double_mat::Zero(nI, nD);
	double_mat H_DD = double_mat::Zero(nD, nD);
	for(size_t k = 0; k < information_rcv.nnz(); k++)
	{	// note only lower triangle is stored in information_rcv
		// but we must fill both lower and upper triangle for eigen inverse
		size_t r = information_rcv.row()[k];
		size_t c = information_rcv.col()[k];
		double v = information_rcv.val()[k];
		//
		bool in_subset = fixed2subset[r] != n_fixed_;
		in_subset     &= fixed2subset[c] != n_fixed_;
		//
		if( in_subset )
		{	size_t i = fixed2subset[r];
			size_t j = fixed2subset[c];
			//
			size_t i_in_I = subset2independent[i];
			size_t j_in_I = subset2independent[j];
			size_t i_in_D = subset2dependent[i];
			size_t j_in_D = subset2dependent[j];
			//
			// H_II
			if( i_in_I < nI && j_in_I < nI )
			{	assert( i_in_D == nD && j_in_D == nD );
				H_II(i_in_I, j_in_I) = v;
				H_II(j_in_I, i_in_I) = v;
			}
			// H_ID
			if( i_in_I < nI && j_in_D < nD )
			{	assert( i_in_D == nD && j_in_I == nI );
				H_ID(i_in_I, j_in_D) = v;
			}
			if( i_in_D < nD && j_in_I < nI )
			{	assert( i_in_I == nI && j_in_D == nD );
				H_ID(j_in_I, i_in_D) = v;
			}
			// H_DD
			if( i_in_D < nD && j_in_D < nD )
			{	assert( i_in_I == nI && j_in_I == nI );
				H_DD(i_in_D, j_in_D) = v;
				H_DD(j_in_D, i_in_D) = v;
			}
		}
	}
	double_mat info_mat = H_II;
	if( nD > 0 )
	{	double_mat H_ID_C   = double_mat( H_ID * C );
		info_mat += H_ID_C.transpose() + H_ID_C + C.transpose() * H_DD * C;
	}
	//
	// create a sparse_mat_info representation of info_mat
	// in column major order
	CppAD::mixed::sparse_mat_info  info_mat_info;
	for(size_t j = 0; j < nI; j++)
	{	for(size_t i = j; i < nI; i++)
		{	double v = info_mat(i, j);
			if( v != 0.0 )
			{	info_mat_info.row.push_back(i);
				info_mat_info.col.push_back(j);
				info_mat_info.val.push_back(v);
			}
		}
	}
	//
	// LDLT factorization of info_mat
	CPPAD_MIXED_LDLT_CLASS ldlt_info_mat(nI);
	ldlt_info_mat.init( info_mat_info );
	bool ok = ldlt_info_mat.update( info_mat_info );
	if( ! ok )
	{	std::string msg =
			"sample_fixed: Implicit information matrix is singular";
		fatal_error(msg);
	}
	size_t negative;
	ldlt_info_mat.logdet(negative);
	if( negative != 0 )
	{	std::string msg = "Implicit information matrix is positive definite";
		warning(msg);
	}
	//
	// -----------------------------------------------------------------------
	// Simulate the samples
	// -----------------------------------------------------------------------
	for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
	{	d_vector w(nI);
		// simulate a normal with mean zero and variance one
		for(size_t k = 0; k < nI; k++)
			w[k] = gsl_ran_gaussian(CppAD::mixed::get_gsl_rng(), 1.0);
		//
		// set v to cholesky factor of info_mat^{-1} times w
		d_vector v(nI);
		ok = ldlt_info_mat.sim_cov(w, v);
		if( ! ok )
		{	std::string msg = "sample_fixed: implicit information matrix"
				" is not positive definite";
			fatal_error(msg);
		}
		//
		// independent variable values
		double_vec alpha_I(nI);
		for(size_t i = 0; i < nI; i++)
			alpha_I[i] = v[i];
		// dependent variables
		double_vec alpha_D = double_vec( C * alpha_I );
		//
		// store in alpha
		double_vec alpha(n_subset);
		for(size_t k = 0; k < nI; k++)
			alpha[I[k]] = alpha_I[k];
		for(size_t k = 0; k < nD; k++)
			alpha[D[k]] = alpha_D[k];
		//
		// store in sample
		for(size_t j = 0; j < n_fixed_; j++)
		{	if( fixed2subset[j] == n_fixed_ )
				sample[ i_sample * n_fixed_ + j] = fixed_opt[j];
			else
			{	size_t k       = fixed2subset[j];
				double fixed_j = fixed_opt[j] + alpha[k];
				//
				// check if this component went out of bounds
				fixed_j = std::min(fixed_j, fixed_upper[j]);
				fixed_j = std::max(fixed_j, fixed_lower[j]);
				//
				// store this component of the sample
				sample[ i_sample * n_fixed_ + j] = fixed_j;
			}
		}
	}
	// -----------------------------------------------------------------------
	return;
}
// -------------------------------------------------------------------------
// BEGIN PROTOTYPE
void cppad_mixed::sample_fixed(
	CppAD::vector<double>&                 sample               ,
	const d_sparse_rcv&                    information_rcv      ,
	const CppAD::mixed::fixed_solution&    solution             ,
	const CppAD::vector<double>&           fixed_lower          ,
	const CppAD::vector<double>&           fixed_upper          ,
	const CppAD::vector<double>&           random_opt           )
// END PROTOTYPE
{	try
	{	try_sample_fixed(
			sample            ,
			information_rcv   ,
			solution          ,
			fixed_lower       ,
			fixed_upper       ,
			random_opt
		);
	}
	catch(const CppAD::mixed::exception& e)
	{	std::string error_message = e.message("sample_fixed");
		fatal_error(error_message);
		assert(false);
	}
	return;
}
