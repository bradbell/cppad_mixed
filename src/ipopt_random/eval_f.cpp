/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_random_eval_f$$
$spell
	eval
	obj
	Ipopt
	Taylor
	vec
	eval
	Ipopt
$$

$section Compute Value of Objective$$

$head Syntax$$
$icode%ok% = eval_f(%n%, %x%, %new_x%, %obj_value%)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the objective
f(x) is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head obj_val$$
set to the initial value of the objective function f(x).

$head ok$$
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head objective_current_$$
After this call, $code objective_current_$$ will be the
value of the objective corresponding to $icode x$$.
Note that if $icode new_x$$ is false,
this value does not change.

$head mixed_object_.ran_like_fun_$$
After this call, the zero order Taylor coefficients in this function
will corresponding to the value of $icode fixed_vec_$$ and
the random effects in $icode x$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_random::eval_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number&         obj_value )  // out
/* %$$
$end
*/
{	try
	{	try_eval_f(n, x, new_x, obj_value);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_random::eval_f");
		for(size_t j = 0; j < n_random_; j++)
			error_random_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_random::try_eval_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number&         obj_value )  // out
{	assert( size_t(n) == n_random_ );
	assert( mixed_object_.ran_like_fun_.Domain() == n_fixed_ + n_random_ );
	assert( mixed_object_.ran_like_fun_.Range()  == 1 );
	//
	if( new_x )
	{	// random effects as a vector
		d_vector random_vec(n_random_);
		for(size_t j = 0; j < n_random_; j++)
			random_vec[j] = x[j];
		//
		// pack both the fixed and random effects into one vector
		d_vector both_vec(n_fixed_ + n_random_);
		mixed_object_.pack(fixed_vec_, random_vec, both_vec);
		//
		// compute the log-density vector
		d_vector vec = mixed_object_.ran_like_fun_.Forward(0, both_vec);
		assert( vec.size() == 1 );
		if( CppAD::hasnan( vec ) ) throw CppAD::mixed::exception(
			"", "objective has a nan"
		);
		// store for re-use
		objective_current_ = vec[0];
	}
	//
	obj_value = objective_current_;
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
