// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin ipopt_nlp_xam$$
$spell
	CppAD
	Ipopt
$$

$section Ipopt Example: Declare Non-linear Program Problem Class$$

$nospell
$srccode%cpp% */
# include <cmath>
# include <coin-or/IpIpoptApplication.hpp>
# include <coin-or/IpTNLP.hpp>
# include <cassert>
namespace {
	// Ipopt types used by this file
	typedef Ipopt::Number                      Number;
	typedef Ipopt::Index                       Index;
	typedef Ipopt::TNLP::IndexStyleEnum        IndexStyleEnum;
	typedef Ipopt::AlgorithmMode               AlgorithmMode;
	typedef Ipopt::IpoptData                   IpoptData;
	typedef Ipopt::IpoptCalculatedQuantities   IpoptCalculatedQuantities;
	//
	// ipopt_nlp_xam
	class ipopt_nlp_xam : public Ipopt::TNLP
	{
	public:
		// factor in the objective function
		const double beta_;
		//
		// did check of solution pass
		bool finalize_solution_ok_;
		//
		// final solution
		std::vector<double> final_solution_;
		//
		// default constructor
		ipopt_nlp_xam(double beta);
		//
		// default destructor
		virtual ~ipopt_nlp_xam(void);
		//
		virtual bool get_nlp_info(
			Index&          n            ,
			Index&          m            ,
			Index&          nnz_jac_g    ,
			Index&          nnz_h_lag    ,
			IndexStyleEnum& index_style
		);
		virtual bool get_bounds_info(
				Index       n        ,
				Number*     x_l      ,
				Number*     x_u      ,
				Index       m        ,
				Number*     g_l      ,
				Number*     g_u
		);
		virtual bool get_starting_point(
			Index           n            ,
			bool            init_x       ,
			Number*         x            ,
			bool            init_z       ,
			Number*         z_L          ,
			Number*         z_U          ,
			Index           m            ,
			bool            init_lambda  ,
			Number*         lambda
		);
		virtual bool eval_f(
			Index           n        ,
			const Number*   x        ,
			bool            new_x    ,
			Number&         obj_value
		);
		virtual bool eval_grad_f(
			Index           n        ,
			const Number*   x        ,
			bool            new_x    ,
			Number*         grad_f
		);
		virtual bool eval_g(
			Index           n        ,
			const Number*   x        ,
			bool            new_x    ,
			Index           m        ,
			Number*         g
		);
		virtual bool eval_jac_g(
			Index           n        ,
			const Number*   x        ,
			bool            new_x    ,
			Index           m        ,
			Index           nele_jac ,
			Index*          iRow     ,
			Index*          jCol     ,
			Number*         values
		);
		virtual bool eval_h(
			Index         n              ,
			const Number* x              ,
			bool          new_x          ,
			Number        obj_factor     ,
			Index         m              ,
			const Number* lambda         ,
			bool          new_lambda     ,
			Index         nele_hess      ,
			Index*        iRow           ,
			Index*        jCol           ,
			Number*       values
		);
		virtual void finalize_solution(
			Ipopt::SolverReturn               status    ,
			Index                             n         ,
			const Number*                     x         ,
			const Number*                     z_L       ,
			const Number*                     z_U       ,
			Index                             m         ,
			const Number*                     g         ,
			const Number*                     lambda    ,
			Number                            obj_value ,
			const Ipopt::IpoptData*           ip_data   ,
			Ipopt::IpoptCalculatedQuantities* ip_cq
		);
		virtual bool intermediate_callback(
			AlgorithmMode               mode,
			Index                       iter,
			Number                      obj_value,
			Number                      inf_pr,
			Number                      inf_du,
			Number                      mu,
			Number                      d_norm,
			Number                      regularization_size,
			Number                      alpha_du,
			Number                      alpha_pr,
			Index                       ls_trials,
			const IpoptData*            ip_data,
			IpoptCalculatedQuantities*  ip_cq
		);
	};
}
/* %$$
$$ $comment end of nospell$$
$end
------------------------------------------------------------------------------
$begin ipopt_xam_ctor$$
$spell
	CppAD
	ipopt_nlp_xam
$$

$section Ipopt Example: Constructor and Destructor$$

$srccode%cpp% */
ipopt_nlp_xam::ipopt_nlp_xam(double beta) : beta_(beta)
{ }
ipopt_nlp_xam::~ipopt_nlp_xam(void)
{ }
/* %$$
$end
------------------------------------------------------------------------------
$begin ipopt_xam_get_nlp_info$$
$spell
	CppAD
	ipopt_xam_get_nlp_info
	nnz_jac
	Jacobian
	bool
	Enum
	bool
	nlp
$$

$section Return Information About Problem Sizes$$

$head Syntax$$
$icode%ok% = get_nlp_info(%n%, %m%, %nnz_jac_g%, %nnz_h_lag%, %index_style%)%$$

$head n$$
is set to the number of variables in the problem (dimension of x).

$head m$$
is set to the number of constraints in the problem (dimension of g(x)).

$head nnz_jac_g$$
is set to the number of nonzero entries in the Jacobian of g(x).

$head nnz_h_lag$$
is set to the number of nonzero entries in the
lower-left triangle of the Hessian of the Lagrangian
$latex f(x) + \lambda^\R{T} g(x)$$.
(The Hessian is symmetric, and hence determined by its lower-left triangle.)

$head index_style$$
is set to the numbering style used for row/col entries in the sparse matrix
format (C_STYLE: 0-based, FORTRAN_STYLE: 1-based).

$head ok$$
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::get_nlp_info(
	Index&          n            ,  // out
	Index&          m            ,  // out
	Index&          nnz_jac_g    ,  // out
	Index&          nnz_h_lag    ,  // out
	IndexStyleEnum& index_style  )  // out
{
	n           = 2;
	m           = 1;
	nnz_jac_g   = 2;
	nnz_h_lag   = 2;
	index_style = C_STYLE;
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_get_bounds_info$$
$spell
	CppAD
	ipopt_nlp_xam
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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::get_bounds_info(
		Index       n        ,   // in
		Number*     x_l      ,   // out
		Number*     x_u      ,   // out
		Index       m        ,   // in
		Number*     g_l      ,   // out
		Number*     g_u      )   // out
{
	assert( n == 2 );
	//
	x_l[0] = 0.0;
	x_u[0] = 2.0;
	//
	x_l[1] = 0.0;
	x_u[1] = 3.0;
	//
	assert( m == 1 );
	//
	g_l[0] = 0.0;
	g_u[0] = 0.0;
	//
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_get_starting_point$$
$spell
	CppAD
	init
	ipopt_nlp_xam
	bool
$$

$section Return Initial Values Where Optimization is Started$$

$head Syntax$$
$icode%ok% = get_starting_point(
	%n%, %init_x%, %x%, %init_z%, %z_L%, %z_U%, %m%, %init_lambda%, %lambda%
)%$$

$head n$$
is the number of variables in the problem (dimension of x).

$head init_x$$
if true, the ipopt options specify that the this routine
will provide an initial value for $icode x$$.

$head x$$
if $icode init_x$$ is true,
set to the initial value for the primal variables (has size $icode n$$).

$head init_z$$
if true, the ipopt options specify that the this routine
will provide an initial value for $icode x$$ upper and lower bound
multipliers.

$head z_L$$
if $icode init_z$$ is true,
set to the initial value for the lower bound multipliers (has size $icode n$$).

$head z_U$$
if $icode init_z$$ is true,
set to the initial value for the upper bound multipliers (has size $icode n$$).

$head init_lambda$$
if true, the ipopt options specify that the this routine
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

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::get_starting_point(
	Index           n            ,  // in
	bool            init_x       ,  // in
	Number*         x            ,  // out
	bool            init_z       ,  // in
	Number*         z_L          ,  // out
	Number*         z_U          ,  // out
	Index           m            ,  // in
	bool            init_lambda  ,  // in
	Number*         lambda       )  // out
{
	assert( n == 2 );
	assert( init_x == true );
	x[0] = 1.0;
	x[1] = 1.0;
	assert( init_z == false );
	assert( m == 1 );
	assert( init_lambda == false );

	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_eval_f$$
$spell
	CppAD
	ipopt_nlp_xam
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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::eval_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number&         obj_value )  // out
{
	assert( n == 2 );
	obj_value  = beta_ * (x[0] - 2.0) * (x[0] - 2.0);
	obj_value += beta_ * (x[1] - 3.0) * (x[1] - 3.0);

	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_eval_grad_f$$
$spell
	CppAD
	ipopt_nlp_xam
	bool
	eval
	const
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

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::eval_grad_f(
	Index           n         ,  // in
	const Number*   x         ,  // in
	bool            new_x     ,  // in
	Number*         grad_f    )  // out
{
	assert( n == 2 );
	grad_f[0] = 2.0 * beta_ * (x[0] - 2.0);
	grad_f[1] = 2.0 * beta_ * (x[1] - 3.0);
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_eval_g$$
$spell
	CppAD
	ipopt_nlp_xam
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
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::eval_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Number*         g        )  // out
{
	assert( n == 2 );
	assert( m == 1 );
	//
	g[0] = x[0] + x[1] - 2.0;
	//
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_eval_jac_g$$
$spell
	CppAD
	ipopt_nlp_xam
	bool
	eval
	const
	nele_jac
	Jacobian
	nnz
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
$cref/nnz_jac_g/ipopt_xam_get_nlp_info/nnz_jac_g/$$.

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

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::eval_jac_g(
	Index           n        ,  // in
	const Number*   x        ,  // in
	bool            new_x    ,  // in
	Index           m        ,  // in
	Index           nele_jac ,  // in
	Index*          iRow     ,  // out
	Index*          jCol     ,  // out
	Number*         values   )  // out
{
	assert( nele_jac == 2 );
	if( values == NULL )
	{
		iRow[0] = 0;
		jCol[0] = 0;
		//
		iRow[1] = 0;
		jCol[1] = 1;
		//
		return true;
	}
	assert( n == 2 );
	assert( m == 1 );
	//
	values[0] = 1.0;
	values[1] = 1.0;
	//
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_eval_h$$
$spell
	CppAD
	ipopt_nlp_xam
	bool
	eval
	const
	obj
	nele_hess
	nnz
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
$cref/nnz_h_lag/ipopt_xam_get_nlp_info/nnz_h_lag/$$.

$head iRow$$
If $icode values$$ is $code NULL$$,
$icode iRow$$ has size $icode nele_hess$$ and is set to the
row indices for the non-zero entries in the
lower (or upper) triangle of the Hessian $latex L_{x,x} (x)$$.

$head jCol$$
If $icode values$$ is $code NULL$$,
$icode jCol$$ has size $icode nele_hess$$ and is set to the
column indices for the non-zero entries in the
lower (or upper) triangle of the Hessian $latex L_{x,x} (x)$$.

$head values$$
If $icode values$$ is not $code NULL$$,
it has size $icode nele_hess$$ and $icode%values%[%k%]%$$
is set to the value of element of the Hessian $latex L_{x,x} (x)$$
with row index $icode%iRow%[%k%]%$$
and column index $icode%jCol%[%k%]%$$.

$head ok$$
if set to false, the optimization will treat this point like
it was not feasible
(the objective function could not be evaluated at this point).

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::eval_h(
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
{
	assert( nele_hess == 2 );
	if( values == NULL )
	{
		iRow[0] = 0;
		jCol[0] = 0;
		//
		iRow[1] = 1;
		jCol[1] = 1;
		//
		return true;
	}
	assert( n == 2 );
	assert( m == 1 );
	//
	values[0] = 2.0 * beta_ * obj_factor;
	values[1] = 2.0 * beta_ * obj_factor;
	//
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_finalize_solution$$
$spell
	CppAD
	ipopt_nlp_xam
	bool
	eval
	const
	obj
	ip
	cq
	namespace
	infeasibility
	doesn't
	Inf
	naninf
	std
	cout
	endl
	fabs
	tol
$$

$section Get Solution Results$$

$head Syntax$$
$codei%finalize_solution(
	%status%, %n%, %x%, %z_L%, %z_U%, %m%, %g%,%$$
$icode%lambda%, %obj_value%, %ip_data%, %ip_cq%
)%$$

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

$head Source$$
$srccode%cpp% */
void ipopt_nlp_xam::finalize_solution(
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
{	bool ok = true;
	using std::fabs;

	// default tolerance
	double tol = 1e-08;

	// check problem dimensions
	ok &= n == 2;
	ok &= m == 1;

	// check that x is feasible
	ok &= (0.0 <= x[0]) && (x[0] <= +2.0);
	ok &= (0.0 <= x[1]) && (x[1] <= +3.0);

	// check that the bound multipliers are feasible
	ok &= (0.0 <= z_L[0]) && (0.0 <= z_L[1]);
	ok &= (0.0 <= z_U[0]) && (0.0 <= z_U[1]);

	// check that the constraint on g(x) is satisfied
	ok &= fabs( x[0] + x[1] - 2.0 ) <= 10. * tol;

	// Check the partial of the Lagrangian w.r.t x[0]
	ok &= fabs(
		2.0 * beta_ * (x[0] - 2.0) + lambda[0] - z_L[0] + z_U[0]
	) <= 10. * tol;

	// Check the partial of the Lagrangian w.r.t x[1]
	ok &= fabs(
		2.0 * beta_ * (x[1] - 3.0) + lambda[0] - z_L[1] + z_U[1]
	) <= 10. * tol;

	// set member variable finalize_solution_ok_
	finalize_solution_ok_ = ok;

	// set member variable final_solution_
	final_solution_.resize(n);
	for(Index j = 0; j < n; j++)
		final_solution_[j] = x[j];
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_xam_intermediate_callback$$
$spell
	Ipopt
	iter
	obj
	inf
	pr
	du
	ls
	ip_cq
	lg
	rg
	optimizer
	mu
	enum
$$

$section Ipopt Example: Optimization Progress Report$$

$head Syntax$$
$icode%ok% = intermediate_callback(
	%mode%,
	%iter%,
	%obj_value%,
	%inf_pr%,
	%inf_du%,
	%mu%,
	%d_norm%,
	%regularization_size%,
	%alpha_du%,
	%alpha_pr%,
	%ls_trials%,
	%ip_data%,
	%ip_cq%
)%$$

$head mode$$
This is an $code enum$$ value and equal to
$code RegularMode$$ or $code RestorationPhaseMode$$.

$head iter$$
See ipopt trace $cref/iter/ipopt_trace/iter/$$.

$head obj_value$$
See ipopt trace $cref/objective/ipopt_trace/objective/$$.

$head inf_pr$$
See ipopt trace $cref/inf_pr/ipopt_trace/inf_pr/$$.

$head inf_du$$
See ipopt trace $cref/inf_du/ipopt_trace/inf_du/$$.

$head mu$$
See ipopt trace $cref/lg(mu)/ipopt_trace/lg(mu)/$$.

$head d_norm$$
See ipopt trace $cref/||d||/ipopt_trace/||d||/$$.

$head regularization_size$$
See ipopt trace $cref/lg(rg)/ipopt_trace/lg(rg)/$$.

$head alpha_du$$
See ipopt trace $cref/alpha_du/ipopt_trace/alpha_du/$$.

$head alpha_pr$$
See ipopt trace $cref/alpha_pr/ipopt_trace/alpha_pr/$$.

$head ls_trials$$
See ipopt trace $cref/ls/ipopt_trace/ls/$$.

$head ok$$
If set to false, the optimizer will terminate with
$code USER_REQUESTED_STOP$$ as the finalize_solution
$cref/status/ipopt_xam_finalize_solution/status/$$.

$head Source$$
$srccode%cpp% */
bool ipopt_nlp_xam::intermediate_callback(
	AlgorithmMode               mode                 ,   // in
	Index                       iter                 ,   // in
	Number                      obj_value            ,   // in
	Number                      inf_pr               ,   // in
	Number                      inf_du               ,   // in
	Number                      mu                   ,   // in
	Number                      d_norm               ,   // in
	Number                      regularization_size  ,   // in
	Number                      alpha_du             ,   // in
	Number                      alpha_pr             ,   // in
	Index                       ls_trials            ,   // in
	const IpoptData*            ip_data              ,   // in
	IpoptCalculatedQuantities*  ip_cq                )   // in
{
	return true;
}
/* %$$
$end
-------------------------------------------------------------------------------
$begin ipopt_run_xam$$
$spell
	CppAD
	ipopt_run_xam
	bool
	eval
	const
	Ptr
	nlp
	sb
	Status status
$$

$section Ipopt: Example and Test$$

$head Syntax$$
$icode%ok% = ipopt_run_xam()%$$

$head ok$$
This return value is true, if the test passes,
and false otherwise.

$head Source$$
$srccode%cpp% */
bool ipopt_run_xam(void)
{	bool ok    = true;
	double tol = 1e-8;
	using Ipopt::SmartPtr;

	// Factor in definition of objective function
	double beta = 3.0;

	// Create an instance of the example problem
	SmartPtr<ipopt_nlp_xam> xam_nlp = new ipopt_nlp_xam(beta);

	// Create an instance of an IpoptApplication
	SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

	// Turn off all Ipopt printed output
	app->Options()->SetIntegerValue("print_level", 0);
	app->Options()->SetStringValue("sb", "yes");
	app->Options()->SetStringValue("derivative_test", "second-order");

	// variable to hold status values returned by app
	Ipopt::ApplicationReturnStatus status;

	// initialize app
	status = app->Initialize();
	ok    &= status == Ipopt::Solve_Succeeded;

	// solve the problem
	status = app->OptimizeTNLP(xam_nlp);
	ok    &= status == Ipopt::Solve_Succeeded;
	ok    &= xam_nlp->finalize_solution_ok_;

	// check the solution
	const std::vector<double>& x(xam_nlp->final_solution_);
	ok    &= x.size() == 2;
	ok    &= fabs( x[0] - 0.5 ) <= 10. * tol;
	ok    &= fabs( x[1] - 1.5 ) <= 10. * tol;

	return ok;
}
/* %$$
$end
*/
