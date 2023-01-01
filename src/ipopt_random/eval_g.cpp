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
$begin ipopt_random_eval_g$$
$spell
	eval
	Ipopt
	Taylor
	vec
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head mixed_object_.ran_like_fun_$$
if $icode new_x$$ is true,
after this call, the zero order Taylor coefficients in this function
will corresponding to the value of $icode fixed_vec_$$ and
the random effects in $icode x$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_random::eval_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Number*         g        )  // out
/* %$$
$end
*/
{	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	//
	if( new_x )
	{	// set the zero order Taylor coefficients in
		// mixed_object_.ran_like_fun_
		Number obj_value;
		eval_f(n, x, new_x, obj_value);
	}
	//
	return true;
}
} } // END_CPPAD_MIXED_NAMESPACE
