// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/configure.hpp>

/*
$begin init_ran_hes$$
$spell
	rcv
	nnz
	CppAD
	init
	cppad
	hes hes
	vec
	const
	Cpp
	logdet
	Cholesky
	namespace
	hpp
	Simplicial
	triangular
	chol
	dismod
	bool
	jac
$$

$section Initialize Hessian of Random Likelihood w.r.t Random Effects$$

$head Syntax$$
$icode%mixed_object%.init_ran_hes(
	%fixed_vec%, %random_vec%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Assumptions$$
The member variables
$cref/init_ran_like_done_/init_ran_like/init_ran_like_done_/$$ and
$cref/init_ran_jac_done_/init_ran_jac/init_ran_jac_done_/$$ are true.

$head init_ran_hes_done_$$
The input value of this member variable must be false.
Upon return it is true.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which the initialization is done.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which the initialization is done.

$head ran_hes_rcv_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_rcv ran_hes_rcv_
%$$
does not matter.
Upon return it contains the
$cref/sparse_rcv/typedef/Sparse Types/sparse_rcv/$$ information
for the lower triangle of the Hessian
$latex \[
	f_{u,u} ( \theta , u )
\]$$
see $cref/f(theta, u)/
	theory/
	Random Likelihood, f(theta, u)
/$$
The matrix is symmetric and hence can be recovered from
its lower triangle.


$head ran_like_fun_$$
If the return value from $cref ran_likelihood_hes$$ is non-empty,
and $code NDEBUG$$ is not defined,
$code ran_likelihood_hes$$ is checked for correctness.

$head Random Effects Index$$
The indices in $code ran_hes_rcv_$$
are relative to both the fixed and random effects with the fixed effects
coming first.
To get the indices relative to just the random effects, subtract
$code n_fixed_$$; e.g, the value
$codei%
	ran_hes_rcv_.val()[%k%]
%$$
corresponds the second partial with respect to the following
two random effect indices:
$codei%
	ran_hes_rcv_.row()[%k%] - n_fixed_
	ran_hes_rcv_.col()[%k%] - n_fixed_
%$$

$head Order$$
The results are in column major order; i.e.,
$codei%
	ran_hes_rcv_.col()[%k%] <= ran_hes_rcv_.col()[%k+1%]
	if( ran_hes_rcv_.col()[%k%] == ran_hes_rcv_.col()[%k+1%] )
		ran_hes_rcv_.row()[%k%] < ran_hes_rcv_.row()[%k+1%]
%$$

$head ran_hes_fun_$$
The input value of the member variables
$codei%
	CppAD::ADFun<double> ran_hes_fun_
%$$
does not matter.
Upon return its zero order forward mode computes
the lower triangle of the sparse Hessian
$latex \[
	f_{u,u} ( \theta , u )
\]$$
in the same order as the elements of
$code ran_hes_rcv_$$ (and $code a1_hes_rcv_$$).

$contents%example/private/ran_hes_fun.cpp
%$$

$end
*/

void cppad_mixed::init_ran_hes(
	const d_vector& fixed_vec     ,
	const d_vector& random_vec    )
{	assert( ! init_ran_hes_done_ );
	assert( init_ran_like_done_ );
	assert( init_ran_jac_done_ );
	//
	size_t n_both = n_fixed_ + n_random_;
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( ran_like_fun_.Domain() == n_both );
	assert( ran_like_fun_.Range()  == 1 );
	//
	// a1_both = (fixed_vec, random_vec)
	a1_vector a1_both(n_both);
	pack(fixed_vec, random_vec, a1_both);
	//
	// a1_w = [ 1.0 ]
	a1_vector a1_w(1);
	a1_w[0]  = 1.0;
	//
	// hes pattern relative to both fixed and random effects
	// (also count number of entries in lower traingle)
	size_t nnz   = ran_jac2hes_rc_.nnz();
	size_t n_low = 0;
	sparse_rc hes_pattern(n_both, n_both, nnz);
	for(size_t k = 0; k < nnz; ++k)
	{	size_t r = ran_jac2hes_rc_.row()[k];
		size_t c = ran_jac2hes_rc_.col()[k];
		assert( r < n_random_ );
		r = r + n_fixed_;
		hes_pattern.set(k, r, c);
		if( r >= c )
			++n_low;
	}
	//
	// subset of sparstiy pattern that we are calculating
	// in column major order
	sparse_rc hes_low(n_random_, n_both, n_low);
	sparse_rc hes_lower(n_both,  n_both, n_low);
	s_vector col_major = hes_pattern.col_major();
	size_t k_low = 0;
	for(size_t k = 0; k < nnz; k++)
	{	size_t ell = col_major[k];
		size_t r   = hes_pattern.row()[ell];
		size_t c   = hes_pattern.col()[ell];
		assert( r >= n_fixed_ );
		if( r >= c )
		{	hes_lower.set(k_low, r, c);
			hes_low.set(k_low, r - n_fixed_, c);
			++k_low;
		}
	}
	assert( k_low == n_low );
	// -----------------------------------------------------------------------
	// create a d_vector containing (theta, u)
	d_vector both(n_both);
	pack(fixed_vec, random_vec, both);
	// -----------------------------------------------------------------------
	// structure used for calculating subset with d_vector results
	ran_hes_rcv_  = sparse_rcv( hes_lower );
	//
	// structure used for calculating subset with a1d_vector results
	a1_sparse_rcv a1_ran_hes_rcv = a1_sparse_rcv( hes_lower );
	//
	// -----------------------------------------------------------------------
	// Declare the independent and dependent variables for taping calculation
	// of Hessian of the random likelihood w.r.t. the random effects
	CppAD::Independent(a1_both);

	// unpack the independent variables into a1_fixed and a1_random
	a1_vector a1_fixed(n_fixed_), a1_random(n_random_);
	unpack(a1_fixed, a1_random, a1_both);

	// Evaluate ran_likelihood_hes. Its row indices are relative
	// to just the random effects, not both fixed and random effects.
	assert( n_low == a1_ran_hes_rcv.nnz() );
	CppAD::vector<size_t> row(n_low), col(n_low);
	for(size_t k = 0; k < n_low; k++)
	{	assert( a1_ran_hes_rcv.row()[k] >= n_fixed_ );
		assert( a1_ran_hes_rcv.col()[k] >= n_fixed_ );
		row[k] = a1_ran_hes_rcv.row()[k] - n_fixed_;
		col[k] = a1_ran_hes_rcv.col()[k] - n_fixed_;
	}
	a1_vector a1_val_out = ran_likelihood_hes(
		a1_fixed, a1_random, row, col
	);
	if( a1_val_out.size() == 0 )
	{	// The user has not defined ran_likelihood_hes, so use AD to calcuate
		// the Hessian of the random likelihood w.r.t the random effects.
		a1_val_out.resize(n_low);
		a1_sparse_rcv a1_subset( hes_low );
		ran_jac_a1fun_.subgraph_jac_rev(a1_both, a1_subset);
		ran_jac_a1fun_.clear_subgraph();
		a1_val_out = a1_subset.val();
	}
# ifndef NDEBUG
	else	check_user_ran_hes(fixed_vec, random_vec);
# endif
	ran_hes_fun_.Dependent(a1_both, a1_val_out);
	ran_hes_fun_.check_for_nan(false);

	// optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	ran_hes_fun_.optimize("no_conditional_skip");
# endif
	//
	init_ran_hes_done_ = true;
}

/*
$begin check_user_ran_hes$$
$spell
	rcv
	cppad
	const
	CppAD
	std
	vec
	hes
	init
$$

$section Check User Defined ran_likelihood_hes$$

$head Syntax$$
$icode%mixed_object%.check_user_ran_hes(%fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.

$head ran_like_fun_$$
The member variable
$codei%
	CppAD::ADFun<double> ran_like_fun_
%$$
must contain a recording of the random likelihood function; see
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$.

$head ran_hes_rcv_$$
The member variable
$codei%
	CppAD::mixed::sparse_rcv  ran_hes_rcv_
%$$
must contain the information for evaluating the random effects
Hessian using $code ran_like_fun_$$; see
$cref/ran_hes_rcv_/init_ran_hes/ran_hes_rcv_/$$.

$head ran_likelihood_hes$$
The return value from $cref ran_likelihood_hes$$ must be non-empty
and is checked for correctness.
To be specific, the return vector is checked using
$cref/ran_like_fun_/Private/ran_like_fun_/$$.
If the results do no agree,
$cref/fatal_error/Public/User Defined Functions/fatal_error/$$ is called with
an appropriate error message.

$end
-----------------------------------------------------------------------------
*/
void cppad_mixed::check_user_ran_hes(
	const d_vector&        fixed_vec       ,
	const d_vector&        random_vec      )
{
	// a1_double versions of fixed_vec and random_vec
	a1_vector a1_fixed(n_fixed_), a1_random(n_random_);
	for(size_t i = 0; i < n_fixed_; i++)
		a1_fixed[i] = fixed_vec[i];
	for(size_t i = 0; i < n_random_; i++)
		a1_random[i] = random_vec[i];

	// row indices in ran_likelihood_hes are relative
	// to just the random effects, not both fixed and random effects.
	size_t K = ran_hes_rcv_.nnz();
	CppAD::vector<size_t> row(K), col(K);
	sparse_rc pattern(n_random_, n_fixed_ + n_random_, K);
	for(size_t k = 0; k < K; k++)
	{	assert( ran_hes_rcv_.row()[k] >= n_fixed_ );
		assert( ran_hes_rcv_.col()[k] >= n_fixed_ );
		row[k] = ran_hes_rcv_.row()[k] - n_fixed_;
		col[k] = ran_hes_rcv_.col()[k] - n_fixed_;
		//
		pattern.set(k, row[k], col[k] + n_fixed_);
	}
	// user's version of Hessian
	a1_vector a1_val_out = ran_likelihood_hes(
		a1_fixed, a1_random, row, col
	);
	assert( a1_val_out.size() != 0 );
	//
	// CppAD's version of Hessian
	a1_sparse_rcv a1_subset(pattern);
	a1_vector a1_both( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, a1_both);
	ran_jac_a1fun_.subgraph_jac_rev(a1_both, a1_subset);
	ran_jac_a1fun_.clear_subgraph();
	//
	double eps = 100. * std::numeric_limits<double>::epsilon();
	bool ok    = a1_val_out.size() == K;
	if( ! ok )
	{	const std::string error_message = "init_ran_hes: "
		"ran_likelihood_hes return value does not have proper size.";
		fatal_error(error_message);
	}
	for(size_t k = 0; k < K; k++)
	{	ok &= CppAD::NearEqual(
			Value(a1_val_out[k]), a1_subset.val()[k], eps, eps
		);
	}
	if( ! ok )
	{	const std::string error_message = "init_ran_hes: "
		"a ran_likelihood_hes Hessian value does not agree with AD value.";
		fatal_error(error_message);
	}
	return;
}
