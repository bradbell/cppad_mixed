/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_eval_f$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	eval
	obj
	const
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
If set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::eval_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number&         obj_value )  // out
/* %$$
$end
*/
{	if( abort_on_eval_error_ )
	{	try_eval_f(n, x, new_x, obj_value);
		obj_value *= scale_f_;
	}
	else
	{	try
		{	try_eval_f(n, x, new_x, obj_value);
			obj_value *= scale_f_;
		}
		catch(const CppAD::mixed::exception& e)
		{	error_message_ = e.message("ipopt_fixed::eval_f");
			for(size_t j = 0; j < n_fixed_; j++)
				error_fixed_[j] = x[j];
			return false;
		}
	}
	return true;
}
void ipopt_fixed::try_eval_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number&         obj_value )  // out
{
	assert( n > 0 && size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	//
	// value of fixed effects corresponding to this x
	for(size_t j = 0; j < n_fixed_; j++)
		fixed_tmp_[j] = double( x[j] );
	//
	// random part of objective
	double H = Number( 0.0 );
	if( n_random_ > 0 )
	{	//
		// compute the optimal random effects corresponding to fixed effects
		if( new_x )
			new_random(fixed_tmp_);
		H = mixed_object_.ran_obj_eval(fixed_tmp_, random_cur_);
	}
	obj_value = Number(H);
	if( fix_likelihood_vec_tmp_.size() == 0 )
		assert( mixed_object_.fix_like_eval(fixed_tmp_).size() == 0 );
	else
	{
		// fixed part of objective
		// (2DO: cache fix_likelihood_vec_tmp_ for eval_g with same x)
		fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_tmp_);
		//
		// only include smooth part of prior in objective
		obj_value += Number( fix_likelihood_vec_tmp_[0] );
		//
		// auxillary variable with index j is constrainted to be
		// greater than absolute value of 1+j component of fixed likelihood
		for(size_t j = 0; j < fix_likelihood_nabs_; j++)
			obj_value += x[n_fixed_ + j];
	}
	if( CppAD::isnan(obj_value) ) throw CppAD::mixed::exception(
		"", "objective is nan"
	);
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
