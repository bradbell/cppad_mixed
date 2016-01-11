// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
/*
$begin ran_like_jac$$
$spell
	Jacobian
	jac
	CppAD
	cppad
	vec
	const
	Cpp
	xam
$$

$section Jacobian of Random Likelihood w.r.t. Random Effects$$

$head Syntax$$
$icode%jac% = %mixed_object%.ran_like_jac( %fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Purpose$$
This routine computes the Jacobian of the random likelihood
$cref/f(theta, u)/theory/Random Likelihood, f(theta, u)/$$
with respect to the random effects vector $latex u$$; i.e.
$latex \[
	f_u ( \theta, u )
\] $$

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

$head jac$$
The return value has prototype
$codei%
	CppAD::vector<a1_double>& %jac%
%$$
It contains the Jacobian $latex f_u ( \theta , u )$$.

$children%
	example/private/ran_like_jac_xam.cpp
%$$
$head Example$$
The file $cref ran_like_jac_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


// ----------------------------------------------------------------------------
// ran_like_jac
CppAD::vector<cppad_mixed::a1_double> cppad_mixed::ran_like_jac(
	const a1d_vector&        fixed_vec   ,
	const a1d_vector&        random_vec  )
{	assert( init_ran_like_done_ );

	// number of fixed and random effects
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_vec.size() );

	// create an a1d_vector containing (theta, u)
	a1d_vector both_vec( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, both_vec);

	// zero order forward mode
	ran_like_a1fun_.Forward(0, both_vec);

	// first order reverse f_{theta,u}^{(1) ( theta , u )
	assert( ran_like_a1fun_.Range() == 1);
	a1d_vector a1_w(1), jac_both(n_fixed_ + n_random_);
	a1_w[0] = 1.0;
	jac_both = ran_like_a1fun_.Reverse(1, a1_w);

	// extract u part of the Jacobian
	a1d_vector jac_ran(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		jac_ran[j] = jac_both[n_fixed_ + j];
	//
	return jac_ran;
}



