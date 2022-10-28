// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_random_eval_jac_g$$
$spell
	Jacobian
	eval
	jac
	nele
	Ipopt
	nnz
	Taylor
	vec
$$

$section Compute Jacobian of Constraint Functions$$

$head Syntax$$
$icode%ok% = eval_jac_g(
	%n%, %x%, %new_x%, %m%, %nele_jac%, %iRow%, %jCol%, %values%
)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the Jacobian
of the constraints $latex \nabla g(x)$$ is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head nele_jac$$
is the number of non-zero elements in the Jacobian of $icode g(x)$$; i.e.,
the same as
$cref/nnz_jac_g/ipopt_random_get_nlp_info/nnz_jac_g/$$.

$head iRow$$
If $icode values$$ is $code NULL$$,
$icode iRow$$ has size $icode nele_jac$$ and is set to the
row indices for the non-zero entries in the Jacobian of the constraints
$latex g_x (x)$$.

$head jCol$$
If $icode values$$ is $code NULL$$,
$icode jCol$$ has size $icode nele_jac$$ and is set to the
column indices for the non-zero entries in the Jacobian of the constraints
$latex g_x (x)$$.

$head values$$
If $icode values$$ is not $code NULL$$,
it has size $icode nele_jac$$ and $icode%values%[%k%]%$$
is set to the value of element of the Jacobian $latex g_x (x)$$
with row index $icode%iRow%[%k%]%$$
and column index $icode%jCol%[%k%]%$$.

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
bool ipopt_random::eval_jac_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Index           nele_jac ,  // in
	Index*          iRow     ,  // out
	Index*          jCol     ,  // out
	Number*         values   )  // out
/* %$$
$end
*/
{	try
	{	try_eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
	}
	catch(const std::exception& e)
	{	error_message_ = "ipopt_random::eval_jac_g: std::exception: ";
		for(size_t j = 0; j < n_random_; j++)
			error_random_[j] = x[j];
		return false;
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_random::eval_jac_g");
		for(size_t j = 0; j < n_random_; j++)
			error_random_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_random::try_eval_jac_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Index           nele_jac ,  // in
	Index*          iRow     ,  // out
	Index*          jCol     ,  // out
	Number*         values   )  // out
{	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	assert( nele_jac == 0 );
	if( values == NULL )
	{	assert( ! new_x );
		return;
	}
	//
	if( new_x )
	{	// set the zero order Taylor coefficients in
		// mixed_object_.ran_like_fun_
		Number obj_value;
		eval_f(n, x, new_x, obj_value);
	}
	return;
}
} } // END_CPPAD_MIXED_NAMESPACE
