/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_random_get_bounds_info$$
$spell
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head Prototype$$
$srccode%cpp% */
bool ipopt_random::get_bounds_info(
		Index       n        ,   // in
		Number*     x_l      ,   // out
		Number*     x_u      ,   // out
		Index       m        ,   // in
		Number*     g_l      ,   // out
		Number*     g_u      )   // out
/* %$$
$end
*/
{
	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	//
	for(size_t j = 0; j < n_random_; j++)
	{	// map infinity to crazy value required by ipopt
		if( random_lower_[j] == - std::numeric_limits<double>::infinity() )
			x_l[j] = nlp_lower_bound_inf_;
		else
			x_l[j] = random_lower_[j];
		//
		if( random_upper_[j] == std::numeric_limits<double>::infinity() )
			x_u[j] = nlp_upper_bound_inf_;
		else
			x_u[j] = random_upper_[j];
	}
	//
	return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
