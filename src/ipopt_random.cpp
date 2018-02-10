// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_random.hpp>
# include <cppad/mixed/exception.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE


/* $$
------------------------------------------------------------------------------
$begin ipopt_random_get_starting_point$$
$spell
	init
	ipopt
$$

$section Return Initial Values Where Optimization is Started$$

$head Syntax$$
$icode%ok% = get_starting_point(
	%n%, %init_x%, %x%, %init_z%, %z_L%, %z_U%, %m%, %init_lambda%, %lambda%
)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head init_x$$
assumed true which means the ipopt options specify that the this routine
will provide an initial value for $icode x$$.

$head x$$
if $icode init_x$$ is true,
set to the initial value for the primal variables (has size $icode n$$).

$head init_z$$
assumes $icode init_z$$ is false.
If it were true, the ipopt options specify that the this routine
will provide an initial value for $icode x$$ upper and lower bound
multipliers.

$head z_L$$
if $icode init_z$$ is true,
set to the initial value for the lower bound multipliers (has size $icode n$$).

$head z_U$$
if $icode init_z$$ is true,
set to the initial value for the upper bound multipliers (has size $icode n$$).

$head init_lambda$$
assumes $icode init_lambda$$ is false.
If it were true, the ipopt options specify that the this routine
will provide an initial value for $icode g(x)$$ upper and lower bound
multipliers.

$head lambda$$
if $icode init_lambda$$ is true,
set to the initial value for the $icode g(x)$$ multipliers
(has size $icode m$$).

$head ok$$
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head Prototype$$
$srccode%cpp% */
bool ipopt_random::get_starting_point(
	Index           n            ,  // in
	bool            init_x       ,  // in
	Number*         x            ,  // out
	bool            init_z       ,  // in
	Number*         z_L          ,  // out
	Number*         z_U          ,  // out
	Index           m            ,  // in
	bool            init_lambda  ,  // in
	Number*         lambda       )  // out
/* %$$
$end
*/
{
	assert( init_x == true );
	assert( init_z == false );
	assert( init_lambda == false );
	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	//
	// initial value for random effects
	for(size_t j = 0; j < n_random_; j++)
		x[j] = random_in_[j];

	return true;
}
/*
-------------------------------------------------------------------------------
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
/*
-------------------------------------------------------------------------------
$begin ipopt_random_eval_grad_f$$
$spell
	eval
	Ipopt
	Taylor
	vec
$$

$section Compute Gradient of the Objective$$

$head Syntax$$
$icode%ok% = eval_grad_f(%n%, %x%, %new_x%, %grad_f%)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the gradient
$latex \nabla f(x)$$ is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head grad_f$$
is set to the value for the gradient $latex \nabla f(x)$$
(has size $icode m$$).

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
bool ipopt_random::eval_grad_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number*         grad_f    )  // out
/* %$$
$end
*/
{	try
	{	try_eval_grad_f(n, x, new_x, grad_f);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_random::eval_grad_f");
		for(size_t j = 0; j < n_random_; j++)
			error_random_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_random::try_eval_grad_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number*         grad_f    )  // out
{	assert( size_t(n) == n_random_ );
	assert( mixed_object_.ran_like_fun_.Domain() == n_fixed_ + n_random_ );
	assert( mixed_object_.ran_like_fun_.Range()  == 1 );
	//
	if( new_x )
	{	// set the zero order Taylor coefficients in
		// mixed_object_.ran_like_fun_
		Number obj_value;
		eval_f(n, x, new_x, obj_value);
	}
	// compute the gradient w.r.t fixed and random effects
	d_vector w(1);
	w[0] = 1.0;
	d_vector dw = mixed_object_.ran_like_fun_.Reverse(1, w);
	if( CppAD::hasnan( dw ) ) throw CppAD::mixed::exception(
		"", "gradient has nan"
	);
	//
	// return gradient w.r.t random effects
	for(size_t j = 0; j < n_random_; j++)
		grad_f[j] = dw[ n_fixed_ + j ];
	return;
}
/*
-------------------------------------------------------------------------------
$begin ipopt_random_eval_g$$
$spell
	eval
	Ipopt
	Taylor
	vec
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
/*
-------------------------------------------------------------------------------
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
/*
-------------------------------------------------------------------------------
$begin ipopt_random_eval_h$$
$spell
	eval
	obj
	nele
	hess
	Ipopt
	nnz
	Teylor
	vec
	hes
	Taylor
$$

$section Compute the Hessian of the Lagrangian$$

$head Syntax$$
$icode%ok% = eval_h(
	%n%, %x%, %new_x%,%obj_factor%, %m%, %lambda%, %new_lambda%,%$$
$icode%nele_hess%, %iRow%, %jCol%, %values%
)%$$

$head Lagrangian$$
The Lagrangian is defined to be
$latex \[
	L(x) = \alpha f(x) + \sum_{i=0}^{m-1} \lambda_i g_i (x)
\] $$

$head mixed_object.quasi_fixed_$$
It is assumed that this member variable is false.

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the value for the primal variables at which the
Hessian of the Lagrangian is computed (has size $icode n$$).

$head new_x$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode x$$.

$head obj_factor$$
is the factor $latex \alpha$$ that multiplies the objective f(x)
in the definition of the Lagrangian.

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head lambda$$
is the value of the constraint multipliers $latex \lambda$$
at which the Hessian is to be evaluated (has size $icode m$$).

$head new_lambda$$
if true, no Ipopt evaluation method was previous called with the same
value for $icode lambda$$.

$head nele_hess$$
is the number of non-zero elements in the Hessian $latex L_{x,x} (x)$$; i.e.,
the same as
$cref/nnz_h_lag/ipopt_random_get_nlp_info/nnz_h_lag/$$.

$head iRow$$
If $icode values$$ is $code NULL$$,
$icode iRow$$ has size $icode nele_hess$$ and is set to the
row indices for the non-zero entries in the
lower triangle of the Hessian $latex L_{x,x} (x)$$.

$head jCol$$
If $icode values$$ is $code NULL$$,
$icode jCol$$ has size $icode nele_hess$$ and is set to the
column indices for the non-zero entries in the
lower triangle of the Hessian $latex L_{x,x} (x)$$.

$head values$$
If $icode values$$ is not $code NULL$$,
it has size $icode nele_hess$$ and $icode%values%[%k%]%$$
is set to the value of element of the Hessian $latex L_{x,x} (x)$$
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

$head mixed_object_.ran_hes_fun_$$
After this call, the zero order Taylor coefficients in this function
will corresponding to the value of $icode fixed_vec_$$ and
the random effects in $icode x$$.

$head Prototype$$
$srccode%cpp% */
bool ipopt_random::eval_h(
	Index         n              ,  // in
	const Number* x              ,  // in
	bool          new_x          ,  // in
	Number        obj_factor     ,  // in
	Index         m              ,  // in
	const Number* lambda         ,  // in
	bool          new_lambda     ,  // in
	Index         nele_hess      ,  // in
	Index*        iRow           ,  // out
	Index*        jCol           ,  // out
	Number*       values         )  // out
/* %$$
$end
*/
{	try
	{	try_eval_h(
			n,
			x,
			new_x,
			obj_factor,
			m,
			lambda,
			new_lambda,
			nele_hess,
			iRow,
			jCol,
			values
		);
	}
	catch(const CppAD::mixed::exception& e)
	{	error_message_ = e.message("ipopt_random::eval_h");
		for(size_t j = 0; j < n_random_; j++)
			error_random_[j] = x[j];
		return false;
	}
	return true;
}
void ipopt_random::try_eval_h(
	Index         n              ,  // in
	const Number* x              ,  // in
	bool          new_x          ,  // in
	Number        obj_factor     ,  // in
	Index         m              ,  // in
	const Number* lambda         ,  // in
	bool          new_lambda     ,  // in
	Index         nele_hess      ,  // in
	Index*        iRow           ,  // out
	Index*        jCol           ,  // out
	Number*       values         )  // out
{	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	assert( size_t(nele_hess) == nnz_h_lag_ );
	//
	if( new_x )
	{	// set the zero order Taylor coefficients in
		// mixed_object_.ran_like_fun_
		Number obj_value;
		eval_f(n, x, new_x, obj_value);
	}
	//
	const s_vector& row( mixed_object_.ran_hes_rcv_.row() );
	const s_vector& col( mixed_object_.ran_hes_rcv_.col() );
	assert( row.size() == nnz_h_lag_ );
	assert( col.size() == nnz_h_lag_ );
	//
	if( values == NULL )
	{	for(size_t k = 0; k < nnz_h_lag_; k++)
		{	// only returning lower triagle of Hessian of objective
			assert( n_fixed_ <= col[k] );
			assert( col[k] <= row[k] );
			assert( row[k] <= n_fixed_ + n_random_ );
			iRow[k] = static_cast<Index>( row[k] - n_fixed_ );
			jCol[k] = static_cast<Index>( col[k] - n_fixed_ );
		}
		assert( ! new_x );
		return;
	}
	//
	// random effects as a vector
	d_vector random_vec(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		random_vec[j] = x[j];
	//
	// pack both the fixed and random effects into one vector
	d_vector both_vec(n_fixed_ + n_random_);
	mixed_object_.pack(fixed_vec_, random_vec, both_vec);
	//
	if( new_x )
	{	// set zero order Taylor coefficient in ran_like_fun_.
		d_vector vec = mixed_object_.ran_like_fun_.Forward(0, both_vec);
		if( CppAD::hasnan( vec ) ) throw CppAD::mixed::exception(
			"", "Hessian has a nan"
		);
	}
	//
	// computes the Hessian of objecive w.r.t random effects  f_uu (theta, u)
	d_vector val = mixed_object_.ran_hes_fun_.Forward(0, both_vec);
	if( CppAD::hasnan( val ) ) throw CppAD::mixed::exception(
		"", "Hessian has a nan"
	);
	assert( val.size() == nnz_h_lag_ );
	//
	// return the values
	for(size_t k = 0; k < nnz_h_lag_; k++)
		values[k] = obj_factor * static_cast<Number>( val[k] );
	//
	return;
}
/*
-------------------------------------------------------------------------------
$begin ipopt_random_finalize_solution$$
$spell
	obj
	ip
	cq
	CppAD
	solution solution
	Ipopt
	namespace
	infeasibility
	doesn't
	Inf
	naninf
$$

$section Get Solution Results$$

$head Syntax$$
$codei%finalize_solution(
	%status%, %n%, %x%, %z_L%, %z_U%, %m%, %g%,%$$
$icode%lambda%, %obj_value%, %ip_data%, %ip_cq%
)%$$

$head solution_$$
This routine checks the solution values and sets the member variable
$codei%
	CppAD::mixed fixed_solution solution_
%$$.

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the final value (best value found) for the primal variables
(has size $icode n$$).

$head z_L$$
is the final value for the $icode x$$ lower bound constraint multipliers
(has size $icode n$$).

$head z_U$$
is the final value for the $icode x$$ upper bound constraint multipliers
(has size $icode n$$).

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head lambda$$
is the final value for the g(x) constraint multipliers $latex \lambda$$.

$head obj_value$$
is the value of the objective f(x) at the final $icode x$$ value.

$head ip_data$$
Unspecified; i.e., not part of the Ipopt user API.

$head ip_cq$$
Unspecified; i.e., not part of the Ipopt user API.

$head status$$
These status values are in the $code Ipopt$$ namespace; e.g.,
$code SUCCESS$$ is short for $code Ipopt::SUCCESS$$:

$subhead SUCCESS$$
Algorithm terminated successfully at a locally optimal point,
satisfying the convergence tolerances (can be specified by options).

$subhead MAXITER_EXCEEDED$$
Maximum number of iterations exceeded (can be specified by an option).

$subhead CPUTIME_EXCEEDED$$
Maximum number of CPU seconds exceeded (can be specified by an option).

$subhead STOP_AT_TINY_STEP$$
Algorithm proceeds with very little progress.

$subhead STOP_AT_ACCEPTABLE_POINT$$
Algorithm stopped at a point that was converged, not to desired
tolerances, but to acceptable tolerances (see the acceptable-... options).

$subhead LOCAL_INFEASIBILITY$$
Algorithm converged to a point of local infeasibility. Problem may be
infeasible.

$subhead USER_REQUESTED_STOP$$
A user call-back function returned false, i.e.,
the user code requested a premature termination of the optimization.

$subhead DIVERGING_ITERATES$$
It seems that the iterates diverge.

$subhead RESTORATION_FAILURE$$
Restoration phase failed, algorithm doesn't know how to proceed.

$subhead ERROR_IN_STEP_COMPUTATION$$
An unrecoverable error occurred while Ipopt tried to compute
the search direction.

$subhead INVALID_NUMBER_DETECTED$$
Algorithm received an invalid number (such as NaN or Inf) from
the NLP; see also option check_derivatives_for_naninf.

$head Prototype$$
$srccode%cpp% */
void ipopt_random::finalize_solution(
	Ipopt::SolverReturn               status    ,  // in
	Index                             n         ,  // in
	const Number*                     x         ,  // in
	const Number*                     z_L       ,  // in
	const Number*                     z_U       ,  // in
	Index                             m         ,  // in
	const Number*                     g         ,  // in
	const Number*                     lambda    ,  // in
	Number                            obj_value ,  // in
	const Ipopt::IpoptData*           ip_data   ,  // in
	Ipopt::IpoptCalculatedQuantities* ip_cq     )  // in
/* %$$
$end
*/
{	assert( size_t(n) == n_random_ );
	assert( m == 0 );
	//
	assert( random_opt_.size() == 0 );
	random_opt_.resize(n_random_);
	for(size_t j = 0; j < n_random_; j++)
		random_opt_[j] = x[j];
	//
	return;
}
/*
---------------------------------------------------------------------------
*/
} } // END_CPPAD_MIXED_NAMESPACE
