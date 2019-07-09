// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-19 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ran_like_hes.hpp>

/*
$begin init_ran_hes$$
$spell
	uu
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
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

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

$head ran_hes_uu_rcv_$$
The input value of the member variable
$codei%
	CppAD::mixed::d_sparse_rcv ran_hes_uu_rcv_
%$$
does not matter.
Upon return $code ran_hes_uu_rcv_.pat()$$ contains the sparsity pattern
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


$subhead Random Effects Index$$
The indices in $icode ran_hes_uu_rcv_$$ are relative to just
the random effects and hence are all less that $code n_random_$$.

$subhead Order$$
The results are in column major order; i.e.,
$codei%
	ran_hes_uu_rcv_.col()[%k%] <= ran_hes_uu_rcv_.col()[%k+1%]
	if( ran_hes_uu_rcv_.col()[%k%] == ran_hes_uu_rcv_.col()[%k+1%] )
		ran_hes_uu_rcv_.row()[%k%] < ran_hes_uu_rcv_.row()[%k+1%]
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
$code ran_hes_uu_rcv_$$ (and $code a1_hes_rcv_$$).

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
	sparse_rc ran_hes_uu_rc(n_random_,  n_random_, n_low);
	s_vector col_major = hes_pattern.col_major();
	size_t k_low = 0;
	for(size_t k = 0; k < nnz; k++)
	{	size_t ell = col_major[k];
		size_t r   = hes_pattern.row()[ell];
		size_t c   = hes_pattern.col()[ell];
		assert( r >= n_fixed_ );
		assert( c >= n_fixed_ );
		if( r >= c )
		{	ran_hes_uu_rc.set(k_low, r - n_fixed_, c - n_fixed_);
			++k_low;
		}
	}
	assert( k_low == n_low );
	//
	// ran_hes_uu_rcv_
	ran_hes_uu_rcv_ = d_sparse_rcv( ran_hes_uu_rc );
	//
	// Declare the independent and dependent variables for taping calculation
	// of Hessian of the random likelihood w.r.t. the random effects
	size_t abort_op_index = 0;
	bool record_compare   = false;
	CppAD::Independent(a1_both, abort_op_index, record_compare);
	//
	// a1_val
	a1_sparse_rcv a1_ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
		n_fixed_, n_random_, ran_jac_a1fun_, ran_hes_uu_rc, a1_both
	);
	const a1_vector& a1_val( a1_ran_hes_uu_rcv.val() );
	//
	// ran_hes_fun_
	ran_hes_fun_.Dependent(a1_both, a1_val);
	ran_hes_fun_.check_for_nan(true);
	//
	// optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	ran_hes_fun_.optimize("no_conditional_skip");
# endif
	//
	init_ran_hes_done_ = true;
}
