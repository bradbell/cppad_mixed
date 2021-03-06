/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_get_bounds_info$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
$$

$section Return Optimization Bounds$$

$head Syntax$$
$icode%ok% = get_bounds_info(%n%, %x_l%, %x_u%, %m%, %g_l%, %g_u%)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x_l$$
set to the lower bounds for $icode x$$ (has size $icode n$$).

$head x_u$$
set to the upper bounds for $icode x$$ (has size $icode n$$).

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head g_l$$
set to the lower bounds for $icode g(x)$$ (has size $icode m$$).

$head g_u$$
set to the upper bounds for $icode g(x)$$ (has size $icode m$$).

$head ok$$
If set to false, the optimization will terminate with status set to
$cref/USER_REQUESTED_STOP
	/ipopt_fixed_finalize_solution/status/USER_REQUESTED_STOP/$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_fixed::get_bounds_info(
		Index       n        ,   // in
		Number*     x_l      ,   // out
		Number*     x_u      ,   // out
		Index       m        ,   // in
		Number*     g_l      ,   // out
		Number*     g_u      )   // out
/* %$$
$end
*/
{	assert( adaptive_called_ );
	//
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );

	for(size_t j = 0; j < n_fixed_; j++)
	{	// map infinity to crazy value required by ipopt
		if( fixed_lower_[j] == - std::numeric_limits<double>::infinity() )
			x_l[j] = nlp_lower_bound_inf_;
		else
			x_l[j] = fixed_lower_[j];
		//
		if( fixed_upper_[j] == std::numeric_limits<double>::infinity() )
			x_u[j] = nlp_upper_bound_inf_;
		else
			x_u[j] = fixed_upper_[j];
	}
	// auxillary varibles for absolute value terms
	for(size_t j = 0; j < fix_likelihood_nabs_; j++)
	{	x_l[n_fixed_ + j] = nlp_lower_bound_inf_;
		x_u[n_fixed_ + j] = nlp_upper_bound_inf_;
	}
	//
	// constraints for absolute value terms
	for(size_t j = 0; j < 2 * fix_likelihood_nabs_; j++)
	{	g_l[j] = 0.0;
		g_u[j] = nlp_upper_bound_inf_;
	}
	//
	// fixed constraints
	for(size_t j = 0; j < n_fix_con_; j++)
	{	size_t i = 2 * fix_likelihood_nabs_ + j;
		//
		g_l[i] = scale_g_[i] * fix_constraint_lower_[j];
		g_u[i] = scale_g_[i] * fix_constraint_upper_[j];
		if( fix_constraint_lower_[j] == fix_constraint_upper_[j] )
			g_u[i] = g_l[i];
	}
	//
	// random constraints
	for(size_t j = 0; j < n_ran_con_; j++)
	{	g_l[2 * fix_likelihood_nabs_ + n_fix_con_ + j] = 0.0;
		g_u[2 * fix_likelihood_nabs_ + n_fix_con_ + j] = 0.0;
	}
	//
	return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
