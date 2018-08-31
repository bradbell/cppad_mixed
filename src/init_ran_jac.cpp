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
# include <cppad/mixed/exception.hpp>
/*
$begin init_ran_jac$$
$spell
	Jacobian
	jac
	CppAD
	cppad
	vec
	const
	Cpp
	init
	hes
	rc
$$

$section Initialize Jacobian of Random Likelihood w.r.t. Random Effects$$

$head Syntax$$
$icode%mixed_object%.init_ran_jac(%fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Assumptions$$
The member variable
$cref/init_ran_like_done_/init_ran_like/init_ran_like_done_/$$ is true.

$head init_ran_jac_done_$$
The input value of this member variable must be false.
Upon return it is true.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<a1_double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<a1_double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.

$children%example/private/ran_jac_fun.cpp
%$$
$head ran_jac_a1fun_$$
The input value of the member variable
$codei%
	CppAD::ADFun<a1_double> ran_jac_a1fun_
%$$
does not matter.
upon return zero order forward mode for this function computes
the Jacobian of the random likelihood  with respect to the random effects;
i.e.,
$latex \[
	f_u ( \theta , u )
\]$$; see
$cref ran_jac_fun.cpp$$.

$head ran_jac2hes_rc_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_rc ran_jac2hes_rc_
%$$
does not matter.
Upon return it contains the sparsity pattern for
$latex f_{uu} ( \theta , u )$$ (not just lower triangle).
The row indices are in the vector $latex u$$; i.e., just the random effects.
The column indices are in the vector  $latex ( \theta , u )$$; i.e.,
both fixed and random effects.

$comment%example/private/ran_jac_fun_.cpp
%$$
----------------------------------------------------------------------------
$end
*/
void cppad_mixed::init_ran_jac(
	const d_vector& fixed_vec     ,
	const d_vector& random_vec    )
{	assert( ! init_ran_jac_done_ );
	assert( init_ran_like_done_ );
	//
	size_t m      = 1;
	size_t n_both = n_fixed_ + n_random_;
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( ran_like_fun_.Domain() == n_both );
	assert( ran_like_fun_.Range()  == m );
	// ------------------------------------------------------------------------
	// ran_jac_a1fun_
	// ------------------------------------------------------------------------
	//
	// a2_both = (fixed_vec, random_vec)
	a2_vector a2_both(n_both);
	pack(fixed_vec, random_vec, a2_both);
	//
	// Record gradient of random likelihood
	size_t abort_op_index = 0;
	bool record_compare   = false;
	CppAD::Independent(a2_both, abort_op_index, record_compare);
	//
	// a2_w = [ 1.0 ]
	a2_vector a2_w(m);
	a2_w[0]  = 1.0;
	//
	// a2_jac
	ran_like_a2fun_.Forward(0, a2_both);
	a2_vector a2_jac_both = ran_like_a2fun_.Reverse(1, a2_w);
	a2_vector a2_jac(n_random_);
	for(size_t j = 0; j < n_random_; ++j)
		a2_jac[j] = a2_jac_both[n_fixed_ + j];
	//
	// ran_jac_a1fun_
	ran_jac_a1fun_.Dependent(a2_both, a2_jac);
# if CPPAD_MIXED_OPTIMIZE_AD_FUNCTION
	ran_jac_a1fun_.optimize()
# endif
	// ------------------------------------------------------------------------
	// ran_jac2hes_rc_
	// ------------------------------------------------------------------------
	//
	// ran_jac_a1fun_ computes f_u ( theta , u ), but we are only interested
	// in the sparsity patter for f_uu ( theta, u ).
	CppAD::vector<bool> select_domain(n_both), select_range(n_random_);
	for(size_t i = 0; i < n_fixed_; i++)
		select_domain[i] = false;             // false for fixed effects
	for(size_t i = 0; i < n_random_; i++)
	{	select_domain[n_fixed_ + i] = true;   // true for random effects
		select_range[i]             = true;
	}
	sparse_rc hes_pattern;
	bool transpose = false;
	ran_jac_a1fun_.subgraph_sparsity(
		select_domain, select_range, transpose, ran_jac2hes_rc_
	);
	assert(ran_jac2hes_rc_.nr() == n_random_);
	assert(ran_jac2hes_rc_.nc() == n_both);
	//
	// ------------------------------------------------------------------------
	init_ran_jac_done_ = true;
}
