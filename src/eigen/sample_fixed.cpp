// $Id:$
/*
-----------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
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
$cref sample_random$$

$head Prototype$$
$srcthisfile%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1%$$

$head Public$$
This $code cppad_mixed$$ $cref base_class$$ member function is public.

$head Purpose$$
This routine draws independent samples from
the asymptotic posterior distribution for the
fixed effects (given the model and the data).

$head Constant Fixed Effects$$
If the upper and lower limits for a fixed effect are equal,
$codei%
	%fixed_lower%[%j%] == %fixed_upper%[%j%]
%$$
we refer to the $th j$$ fixed effect as constant.

$head Constraints$$
Only the constant fixed effect constants are taken into account.
This can result in samples that are outside the limits
for the fixed effects that are not constant.
You may want to adjust the samples to be within
their upper and lower limits before you use them.

$head Covariance$$
Each sample of the fixed effects
(excluding the constant fixed effects)
has covariance equal to the inverse of the
$cref/information matrix/information_mat/$$
(where the constant fixed effects have been removed).

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
and for fixed $latex i$$, they have the specified
$cref/covariance/sample_fixed/Covariance/$$.

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

$children%example/user/sample_fixed.cpp
	%src/eigen/sample_conditional.cpp
%$$
$head Example$$
The file $cref sample_fixed.cpp$$ is an example
and test of $code sample_fixed$$.

$head Other Method$$
The routine $cref sample_conditional$$, is a old method that is no
longer used for computing these samples.

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
	//
	// solution
	assert( solution.fixed_opt.size() == n_fixed_ );
	assert( solution.fixed_lag.size() == n_fixed_ );
	assert( solution.fix_con_lag.size() == n_fix_con );
	assert( solution.ran_con_lag.size() == A_rcv_.nr() );
	//
	// random_opt
	assert( random_opt.size() == n_random_ );
	//
	// number of samples
	size_t n_sample = sample.size() / n_fixed_;
	//
	// optimal fixed effects
	const d_vector& fixed_opt( solution.fixed_opt );
	//
	// Determine the subset of variables that do not have lower equal to upper
	CppAD::vector<size_t> fixed2subset(n_fixed_);
	size_t n_subset = 0;
	for(size_t j = 0; j < n_fixed_; j++)
	{	if( fixed_lower[j] == fixed_upper[j] )
			fixed2subset[j] = n_fixed_;
		else
			fixed2subset[j] = n_subset++;
	}
	assert( n_subset <= n_fixed_ );
	//
	// create a sparse_rcv representation of information matrix on subset
	// and in column major order
	CppAD::vector<size_t> col_major = information_rcv.col_major();
	size_t count = 0;
	for(size_t ell = 0; ell < information_rcv.nnz(); ++ell)
	{	size_t k = col_major[ell];
		size_t i = information_rcv.row()[k];
		size_t j = information_rcv.col()[k];
		i        = fixed2subset[i];
		j        = fixed2subset[j];
		if( i != n_fixed_ && j != n_fixed_ && j <= i )
		{	assert( i < n_subset && j < n_subset );
			++count;
		}
	}
	sparse_rc info_mat_rc(n_subset, n_subset, count);
	count = 0;
	for(size_t ell = 0; ell < information_rcv.nnz(); ++ell)
	{	size_t k = col_major[ell];
		size_t i = information_rcv.row()[k];
		size_t j = information_rcv.col()[k];
		i        = fixed2subset[i];
		j        = fixed2subset[j];
		if( i != n_fixed_ && j != n_fixed_ && j <= i )
		{	assert( i < n_subset && j < n_subset );
			info_mat_rc.set(count++, i, j);
		}
	}
	assert( count == info_mat_rc.nnz() );
	d_sparse_rcv info_mat_rcv( info_mat_rc );
	count = 0;
	for(size_t ell = 0; ell < information_rcv.nnz(); ++ell)
	{	size_t k = col_major[ell];
		size_t i = information_rcv.row()[k];
		size_t j = information_rcv.col()[k];
		double v = information_rcv.val()[k];
		i        = fixed2subset[i];
		j        = fixed2subset[j];
		if( i != n_fixed_ && j != n_fixed_ && j <= i )
		{	assert( i < n_subset && j < n_subset );
			info_mat_rcv.set(count++, v);
		}
	}
	assert( count == info_mat_rcv.nnz() );
	//
	// LDLT factorization of info_mat
	CPPAD_MIXED_LDLT_CLASS ldlt_info_mat(n_subset);
	ldlt_info_mat.init( info_mat_rcv.pat() );
	bool ok = ldlt_info_mat.update( info_mat_rcv );
	if( ! ok )
	{	std::string msg =
			"sample_fixed: information matrix is singular";
		fatal_error(msg);
	}
	size_t negative;
	ldlt_info_mat.logdet(negative);
	if( negative != 0 )
	{	std::string msg =
			"sample_fixed: information matrix is not positive definite";
		warning(msg);
	}
	//
	// -----------------------------------------------------------------------
	// Simulate the samples
	// -----------------------------------------------------------------------
	for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
	{	d_vector w(n_subset);
		// simulate a normal with mean zero and variance one
		for(size_t k = 0; k < n_subset; k++)
			w[k] = gsl_ran_gaussian(CppAD::mixed::get_gsl_rng(), 1.0);
		//
		// set v to cholesky factor of info_mat^{-1} times w
		d_vector v(n_subset);
		ok = ldlt_info_mat.sim_cov(w, v);
		if( ! ok )
		{	std::string msg = "sample_fixed: implicit information matrix"
				" is not positive definite";
			fatal_error(msg);
		}
		//
		// store in sample
		for(size_t j = 0; j < n_fixed_; j++)
		{	if( fixed2subset[j] == n_fixed_ )
				sample[ i_sample * n_fixed_ + j] = fixed_lower[j];
			else
			{	size_t k       = fixed2subset[j];
				//
				// store this component of the sample
				sample[ i_sample * n_fixed_ + j] = fixed_opt[j] + v[k];
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
