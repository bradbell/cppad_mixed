/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_eval_g$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	const
	eval
$$

$section Compute Value of Constraint Functions$$

$head Syntax$$
$icode%ok% = eval_g(%n%, %x%, %new_x%, %m%, %g%)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the constraints
$latex g(x)$$ is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head g$$
is set to the value for the constraint functions (has size $icode m$$).

$head ok$$
if set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::eval_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Number*         g        )  // out
/* %$$
$end
*/
{	if( abort_on_eval_error_ )
	{	try_eval_g(n, x, new_x, m, g);
		for(size_t i = 0; i < size_t(m); i++)
			g[i] *= scale_g_[i];
	}
	else
	{	try
		{	try_eval_g(n, x, new_x, m, g);
			for(size_t i = 0; i < size_t(m); i++)
				g[i] *= scale_g_[i];
		}
		catch(const std::exception& e)
		{	error_message_ = "ipopt_fixed::eval_g: std::exception: ";
			for(size_t j = 0; j < n_fixed_; j++)
				error_fixed_[j] = x[j];
			return false;
		}
		catch(const CppAD::mixed::exception& e)
		{	error_message_ = e.message("ipopt_fixed::eval_g");
			for(size_t j = 0; j < n_fixed_; j++)
				error_fixed_[j] = x[j];
			return false;
		}
	}
	return true;
}
void ipopt_fixed::try_eval_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Number*         g        )  // out
{
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	//
	// fixed effects
	for(size_t j = 0; j < n_fixed_; j++)
		fixed_tmp_[j] = double( x[j] );
	//
	// check if this is a new x
	if( new_x && n_random_ > 0 )
		new_random(fixed_tmp_);
	//
	// fixed likelihood
	// (2DO: cache fix_likelihood_vec_tmp_ for eval_f with same x)
	fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_tmp_);
	//
	// constraint part of absolute value terms in fixed likelihood
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
	{	// x[n_fixed_ + j] >= fix_likelihood_vec_tmp_[1 + j];
		assert( 2 * j + 1 < size_t(m) );
		// g[2 * j] >= 0
		g[2 * j] = Number(x[n_fixed_ + j] - fix_likelihood_vec_tmp_[1 + j]);
		// x[n_fixed_ + j] >= - fix_likelihood_vec_tmp_[1 + j]
		// g[2 * j + 1] >= 0
		g[2*j+1] = Number(x[n_fixed_ + j] + fix_likelihood_vec_tmp_[1 + j]);
	}
	//
	// fixed constraints
	assert( c_vec_tmp_.size() == n_fix_con_ );
	c_vec_tmp_ = mixed_object_.fix_con_eval(fixed_tmp_);
	for(size_t j = 0; j < n_fix_con_; j++)
	{	assert( 2 * fix_likelihood_nabs_ + j < size_t(m) );
		g[2 * fix_likelihood_nabs_ + j] = c_vec_tmp_[j];
	}
	//
	// random constraints
	assert( A_uhat_tmp_.size() == n_ran_con_ );
	mixed_object_.ran_con_eval(random_cur_, A_uhat_tmp_);
	for(size_t j = 0; j < n_ran_con_; j++)
	{	assert( 2 * fix_likelihood_nabs_ + n_fix_con_ + j < size_t(m) );
		g[2 * fix_likelihood_nabs_ + n_fix_con_ + j] = A_uhat_tmp_[j];
	}
	for(size_t i = 0; i < size_t(m); i++)
	{	if( CppAD::isnan( g[i] ) ) throw CppAD::mixed::exception(
			"", "constaint function has a nan"
		);
	}
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
