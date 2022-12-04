// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
{xrst_begin ipopt_nlp_xam}

Ipopt Example: Declare Non-linear Program Problem Class
#######################################################

{xrst_spell_off}
{xrst_code cpp} */
# include <cmath>
# include <iomanip>
# include <coin-or/IpIpoptApplication.hpp>
# include <coin-or/IpTNLP.hpp>
# include <cassert>
// needed by intermediate_callback to determine if in restoration mode
# include <coin-or/IpIpoptCalculatedQuantities.hpp>
# include <coin-or/IpOrigIpoptNLP.hpp>

# define PRINT_TRACE 0
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
      // final solution: primal variables
      std::vector<double> final_solution_x_;
      // final solution: lower bound Lagrange multipliers
      std::vector<double> final_solution_z_L_;
      // final solution: upper bound Lagrange multipliers
      std::vector<double> final_solution_z_U_;
      // final solution: g(x) constraint Lagrange multipliers
      std::vector<double> final_solution_lambda_;
      //
      // default constructor
      ipopt_nlp_xam(double beta);
      //
      // default destructor
      virtual ~ipopt_nlp_xam(void);
      //
      bool get_nlp_info(
         Index&          n            ,
         Index&          m            ,
         Index&          nnz_jac_g    ,
         Index&          nnz_h_lag    ,
         IndexStyleEnum& index_style
      ) override;
      bool get_bounds_info(
            Index       n        ,
            Number*     x_l      ,
            Number*     x_u      ,
            Index       m        ,
            Number*     g_l      ,
            Number*     g_u
      ) override;
      bool get_starting_point(
         Index           n            ,
         bool            init_x       ,
         Number*         x            ,
         bool            init_z       ,
         Number*         z_L          ,
         Number*         z_U          ,
         Index           m            ,
         bool            init_lambda  ,
         Number*         lambda
      ) override;
      bool eval_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number&         obj_value
      ) override;
      bool eval_grad_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number*         grad_f
      ) override;
      bool eval_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Number*         g
      ) override;
      bool eval_jac_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Index           nele_jac ,
         Index*          iRow     ,
         Index*          jCol     ,
         Number*         values
      ) override;
      bool eval_h(
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
      ) override;
      void finalize_solution(
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
      ) override;
      bool intermediate_callback(
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
      ) override;
   };
}
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_nlp_xam}
------------------------------------------------------------------------------
{xrst_begin ipopt_xam_ctor}

Ipopt Example: Constructor and Destructor
#########################################

{xrst_spell_off}
{xrst_code cpp} */
ipopt_nlp_xam::ipopt_nlp_xam(double beta) : beta_(beta)
{ }
ipopt_nlp_xam::~ipopt_nlp_xam(void)
{ }
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_ctor}
------------------------------------------------------------------------------
{xrst_begin ipopt_xam_get_nlp_info}
{xrst_spell
   nlp
   nnz
   nonzero
}

Return Information About Problem Sizes
######################################

Syntax
******
*ok* = ``get_nlp_info`` ( *n* , *m* , *nnz_jac_g* , *nnz_h_lag* , *index_style* )

n
*
is set to the number of variables in the problem (dimension of x).

m
*
is set to the number of constraints in the problem (dimension of g(x)).

nnz_jac_g
*********
is set to the number of nonzero entries in the Jacobian of g(x).

nnz_h_lag
*********
is set to the number of nonzero entries in the
lower-left triangle of the Hessian of the Lagrangian
:math:`f(x) + \lambda^\R{T} g(x)`.
(The Hessian is symmetric, and hence determined by its lower-left triangle.)

index_style
***********
is set to the numbering style used for row/col entries in the sparse matrix
format (C_STYLE: 0-based, FORTRAN_STYLE: 1-based).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_get_nlp_info}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_get_bounds_info}

Return Optimization Bounds
##########################

Syntax
******
*ok* = ``get_bounds_info`` ( *n* , *x_l* , *x_u* , *m* , *g_l* , *g_u* )

n
*
is the number of variables in the problem (dimension of x).

x_l
***
set to the lower bounds for *x* (has size *n* ).

x_u
***
set to the upper bounds for *x* (has size *n* ).

m
*
is the number of constraints in the problem (dimension of g(x)).

g_l
***
set to the lower bounds for *g* ( *x* ) (has size *m* ).

g_u
***
set to the upper bounds for *g* ( *x* ) (has size *m* ).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_get_bounds_info}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_get_starting_point}

Return Initial Values Where Optimization is Started
###################################################

Syntax
******

| *ok* = ``get_starting_point`` (
| |tab| *n* , *init_x* , *x* , *init_z* , *z_L* , *z_U* , *m* , *init_lambda* , *lambda*
| )

n
*
is the number of variables in the problem (dimension of x).

init_x
******
if true, the ipopt options specify that the this routine
will provide an initial value for *x* .

x
*
if *init_x* is true,
set to the initial value for the primal variables (has size *n* ).

init_z
******
if true, the ipopt options specify that the this routine
will provide an initial value for *x* upper and lower bound
multipliers.
This is true during a warm start and false otherwise.

z_L
***
if *init_z* is true,
set to the initial value for the lower bound multipliers (has size *n* ).

z_U
***
if *init_z* is true,
set to the initial value for the upper bound multipliers (has size *n* ).

init_lambda
***********
if true, the ipopt options specify that the this routine
will provide an initial value for *g* ( *x* ) upper and lower bound
multipliers.
This is true during a warm start and false otherwise.

lambda
******
if *init_lambda* is true,
set to the initial value for the *g* ( *x* ) multipliers
(has size *m* ).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
   assert( m == 1 );
   bool warm_start = init_z || init_lambda;
   if( warm_start )
   {  // use previous final solution as starting point
      assert( init_z == true );
      assert( init_lambda == true );
      assert( final_solution_x_.size()   == size_t(n) );
      assert( final_solution_z_L_.size() == size_t(n) );
      assert( final_solution_z_U_.size() == size_t(n) );
      assert( final_solution_lambda_.size() == size_t(m) );
      for(Index j = 0; j < n; ++j)
      {  x[j]   = final_solution_x_[j];
         z_L[j] = final_solution_z_L_[j];
         z_U[j] = final_solution_z_U_[j];
      }
      for(Index i = 0; i < m; ++i)
         lambda[i] = final_solution_lambda_[i];
   }
   else
   {  assert( init_z == false );
      assert( init_lambda == false );
      x[0] = 1.0;
      x[1] = 1.0;
   }
   return true;
}
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_get_starting_point}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_eval_f}

Compute Value of Objective
##########################

Syntax
******
*ok* = ``eval_f`` ( *n* , *x* , *new_x* , *obj_value* )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the objective
f(x) is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

obj_val
*******
set to the initial value of the objective function f(x).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_eval_f}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_eval_grad_f}

Compute Gradient of the Objective
#################################

Syntax
******
*ok* = ``eval_grad_f`` ( *n* , *x* , *new_x* , *grad_f* )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the gradient
:math:`\nabla f(x)` is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

grad_f
******
is set to the value for the gradient :math:`\nabla f(x)`
(has size *m* ).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_eval_grad_f}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_eval_g}

Compute Value of Constraint Functions
#####################################

Syntax
******
*ok* = ``eval_g`` ( *n* , *x* , *new_x* , *m* , *g* )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the value for the primal variables at which the constraints
:math:`g(x)` is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

m
*
is the number of constraints in the problem (dimension of g(x)).

g
*
is set to the value for the constraint functions (has size *m* ).

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_eval_g}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_eval_jac_g}
{xrst_spell
   nele
}

Compute Jacobian of Constraint Functions
########################################

Syntax
******

| *ok* = ``eval_jac_g`` (
| |tab| *n* , *x* , *new_x* , *m* , *nele_jac* , *iRow* , *jCol* , *values*
| )

n
*
is the number of variables in the problem (dimension of x).

x
*
If *values* is ``NULL`` , *x* should also be ``NULL`` .
Otherwise it is the value for the primal variables at which the Jacobian
of the constraints :math:`\nabla g(x)` is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

m
*
is the number of constraints in the problem (dimension of g(x)).

nele_jac
********
is the number of non-zero elements in the Jacobian of *g* ( *x* ) ; i.e.,
the same as
:ref:`ipopt_xam_get_nlp_info@nnz_jac_g` .

iRow
****
If *values* is ``NULL`` ,
*iRow* has size *nele_jac* and is set to the
row indices for the non-zero entries in the Jacobian of the constraints
:math:`g_x (x)`.

jCol
****
If *values* is ``NULL`` ,
*jCol* has size *nele_jac* and is set to the
column indices for the non-zero entries in the Jacobian of the constraints
:math:`g_x (x)`.

values
******
If *values* is not ``NULL`` ,
it has size *nele_jac* and *values* [ *k* ]
is set to the value of element of the Jacobian :math:`g_x (x)`
with row index *iRow* [ *k* ]
and column index *jCol* [ *k* ] .

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_eval_jac_g}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_eval_h}
{xrst_spell
   hess
   nele
}

Compute the Hessian of the Lagrangian
#####################################

Syntax
******

| *ok* = ``eval_h`` (
| |tab| *n* , *x* , *new_x* , *obj_factor* , *m* , *lambda* , *new_lambda* , *nele_hess* , *iRow* , *jCol* , *values*
| )

Lagrangian
**********
The Lagrangian is defined to be

.. math::

   L(x) = \alpha f(x) + \sum_{i=0}^{m-1} \lambda_i g_i (x)

n
*
is the number of variables in the problem (dimension of x).

x
*
If *values* is ``NULL`` , *x* should also be ``NULL`` .
Otherwise it is the value for the primal variables at which the
Hessian of the Lagrangian is computed (has size *n* ).

new_x
*****
if true, no Ipopt evaluation method was previous called with the same
value for *x* .

obj_factor
**********
is the factor :math:`\alpha` that multiplies the objective f(x)
in the definition of the Lagrangian.

m
*
is the number of constraints in the problem (dimension of g(x)).

lambda
******
is the value of the constraint multipliers :math:`\lambda`
at which the Hessian is to be evaluated (has size *m* ).

new_lambda
**********
if true, no Ipopt evaluation method was previous called with the same
value for *lambda* .

nele_hess
*********
is the number of non-zero elements in the Hessian :math:`L_{x,x} (x)`; i.e.,
the same as
:ref:`ipopt_xam_get_nlp_info@nnz_h_lag` .

iRow
****
If *values* is ``NULL`` ,
*iRow* has size *nele_hess* and is set to the
row indices for the non-zero entries in the
lower (or upper) triangle of the Hessian :math:`L_{x,x} (x)`.

jCol
****
If *values* is ``NULL`` ,
*jCol* has size *nele_hess* and is set to the
column indices for the non-zero entries in the
lower (or upper) triangle of the Hessian :math:`L_{x,x} (x)`.

values
******
If *values* is not ``NULL`` ,
it has size *nele_hess* and *values* [ *k* ]
is set to the value of element of the Hessian :math:`L_{x,x} (x)`
with row index *iRow* [ *k* ]
and column index *jCol* [ *k* ] .

ok
**
if set to false, the optimization will treat this point like
it was not feasible
(the objective function could not be evaluated at this point).

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_eval_h}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_finalize_solution}
{xrst_spell
   cq
   infeasibility
   infeasible
   iterates
   namespace
   naninf
   unrecoverable
}

Get Solution Results
####################

Syntax
******

| ``finalize_solution`` (
| |tab| *status* , *n* , *x* , *z_L* , *z_U* , *m* , *g* , *lambda* , *obj_value* , *ip_data* , *ip_cq*
| )

n
*
is the number of variables in the problem (dimension of x).

x
*
is the final value (best value found) for the primal variables
(has size *n* ).

z_L
***
is the final value for the *x* lower bound constraint multipliers
(has size *n* ).

z_U
***
is the final value for the *x* upper bound constraint multipliers
(has size *n* ).

m
*
is the number of constraints in the problem (dimension of g(x)).

lambda
******
is the final value for the g(x) constraint multipliers :math:`\lambda`.

obj_value
*********
is the value of the objective f(x) at the final *x* value.

ip_data
*******
Unspecified; i.e., not part of the Ipopt user API.

ip_cq
*****
Unspecified; i.e., not part of the Ipopt user API.

status
******
These status values are in the ``Ipopt`` namespace; e.g.,
``SUCCESS`` is short for ``Ipopt::SUCCESS`` :

SUCCESS
=======
Algorithm terminated successfully at a locally optimal point,
satisfying the convergence tolerances (can be specified by options).

MAXITER_EXCEEDED
================
Maximum number of iterations exceeded (can be specified by an option).

CPUTIME_EXCEEDED
================
Maximum number of CPU seconds exceeded (can be specified by an option).

STOP_AT_TINY_STEP
=================
Algorithm proceeds with very little progress.

STOP_AT_ACCEPTABLE_POINT
========================
Algorithm stopped at a point that was converged, not to desired
tolerances, but to acceptable tolerances (see the acceptable-... options).

LOCAL_INFEASIBILITY
===================
Algorithm converged to a point of local infeasibility. Problem may be
infeasible.

USER_REQUESTED_STOP
===================
A user call-back function returned false, i.e.,
the user code requested a premature termination of the optimization.

DIVERGING_ITERATES
==================
It seems that the iterates diverge.

RESTORATION_FAILURE
===================
Restoration phase failed, algorithm doesn't know how to proceed.

ERROR_IN_STEP_COMPUTATION
=========================
An unrecoverable error occurred while Ipopt tried to compute
the search direction.

INVALID_NUMBER_DETECTED
=======================
Algorithm received an invalid number (such as NaN or Inf) from
the NLP; see also option check_derivatives_for_naninf.

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
{  bool ok = true;
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

   // set member variable corresponding to final solution
   final_solution_x_.resize(n);
   final_solution_z_L_.resize(n);
   final_solution_z_U_.resize(n);
   final_solution_lambda_.resize(m);
   final_solution_x_.resize(n);
   for(Index j = 0; j < n; ++j)
   {  final_solution_x_[j]     = x[j];
      final_solution_z_L_[j]   = z_L[j];
      final_solution_z_U_[j]   = z_U[j];
   }
   for(Index i = 0; i < m; ++i)
      final_solution_lambda_[i] = lambda[i];
}
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_finalize_solution}
-------------------------------------------------------------------------------
{xrst_begin ipopt_xam_intermediate_callback}
{xrst_spell
   cq
   enum
   iter
   ls
   optimizer
   pr
}

Ipopt Example: Optimization Progress Report
###########################################

Syntax
******

| *ok* = ``intermediate_callback`` (
| |tab| *mode* ,
| |tab| *iter* ,
| |tab| *obj_value* ,
| |tab| *inf_pr* ,
| |tab| *inf_du* ,
| |tab| *mu* ,
| |tab| *d_norm* ,
| |tab| *regularization_size* ,
| |tab| *alpha_du* ,
| |tab| *alpha_pr* ,
| |tab| *ls_trials* ,
| |tab| *ip_data* ,
| |tab| *ip_cq*
| )

mode
****
This is an ``enum`` value and equal to
``RegularMode`` or ``RestorationPhaseMode`` .

iter
****
See ipopt trace :ref:`ipopt_trace@iter` .

obj_value
*********
See ipopt trace :ref:`ipopt_trace@objective` .

inf_pr
******
See ipopt trace :ref:`ipopt_trace@inf_pr` .

inf_du
******
See ipopt trace :ref:`ipopt_trace@inf_du` .

mu
**
See ipopt trace :ref:`ipopt_trace@lg(mu)` .

d_norm
******
See ipopt trace :ref:`ipopt_trace@||d||` .

regularization_size
*******************
See ipopt trace :ref:`ipopt_trace@lg(rg)` .

alpha_du
********
See ipopt trace :ref:`ipopt_trace@alpha_du` .

alpha_pr
********
See ipopt trace :ref:`ipopt_trace@alpha_pr` .

ls_trials
*********
See ipopt trace :ref:`ipopt_trace@ls` .

ok
**
If set to false, the optimizer will terminate with
``USER_REQUESTED_STOP`` as the finalize_solution
:ref:`ipopt_xam_finalize_solution@status` .

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
# if PRINT_TRACE
   using std::cout;
   using std::right;
   using std::setw;
   //
   bool restoration = Ipopt::GetRawPtr(ip_cq->GetIpoptNLP()) == nullptr;
   //
   std::streamsize restore = cout.precision();
   cout.precision(2);
   cout << std::scientific;
   if( iter % 10 == 0 )
   {  cout << right << setw(4) << "iter";
      cout << right << setw(10) << "objective";
      cout << right << setw(10) << "inf_pr";
      cout << right << setw(10) << "inf_du";
      cout << right << setw(10) << "||d||";
      cout << right << setw(6)  << "l(mu)";
      cout << right << setw(6)  << "l(rg)";
      cout << right << setw(10) << "alpha_du";
      cout << right << setw(10) << "alpha_pr";
      cout << right << setw(3)  << "ls";
      cout << right << setw(3)  << "r";
      cout << "\n";
   }
   cout << right << setw(4) << iter;
   cout << right << setw(10) << obj_value;
   cout << right << setw(10) << inf_pr;
   cout << right << setw(10) << inf_du;
   cout << right << setw(10) << d_norm;
   cout.precision(0);
   cout << std::fixed;
   cout << right << setw(6) << std::log10(mu);
   if( regularization_size <= 0.0 )
      cout << right << setw(6) << "-";
   else
      cout << right << setw(6) << std::log10( regularization_size );
   cout << std::scientific;
   cout.precision(2);
   cout << right << setw(10) << alpha_du;
   cout << right << setw(10) << alpha_pr;
   cout << right << setw(3) << ls_trials;
   cout << right << setw(3) << restoration;
   //
   cout << "\n";
   cout << std::defaultfloat;
   cout.precision(restore);
# endif
   //
   return true;
}
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_xam_intermediate_callback}
-------------------------------------------------------------------------------
{xrst_begin ipopt_run_xam}

Ipopt: Example and Test
#######################

Syntax
******
*ok* = ``ipopt_run_xam`` ()

ok
**
This return value is true, if the test passes,
and false otherwise.

Source
******
{xrst_spell_off}
{xrst_code cpp} */
bool ipopt_run_xam(void)
{  bool ok    = true;
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
   app->Options()->SetIntegerValue("max_iter", 4);
   app->Options()->SetStringValue("sb", "yes");
   app->Options()->SetStringValue("derivative_test", "second-order");

   // variable to hold status values returned by app
   Ipopt::ApplicationReturnStatus status;

   // initialize app
   status = app->Initialize();
   ok    &= status == Ipopt::Solve_Succeeded;

   // first attempt should fail because max_iter is so small
   status = app->OptimizeTNLP(xam_nlp);
   ok    &= status != Ipopt::Solve_Succeeded;
   //
   // second attemp should succeed with a warm start
   app->Options()->SetStringValue("warm_start_init_point", "yes");
   status = app->OptimizeTNLP(xam_nlp);
   //
   ok    &= status == Ipopt::Solve_Succeeded;
   ok    &= xam_nlp->finalize_solution_ok_;

   // check solution primal variables x
   const std::vector<double>& x(xam_nlp->final_solution_x_);
   ok    &= x.size() == 2;
   ok    &= std::fabs( x[0] - 0.5 ) <= 10. * tol;
   ok    &= std::fabs( x[1] - 1.5 ) <= 10. * tol;

   // check solution g(x) Lagrange multipliers lambda
   const std::vector<double>& lambda(xam_nlp->final_solution_lambda_);
   ok    &= lambda.size() == 1;
   ok    &= std::fabs(lambda[0] - 3.0 * beta) <= 10. * tol;

   // check solution bound constraint Lagrange multipliers z_L and z_U
   const std::vector<double>& z_L(xam_nlp->final_solution_z_L_);
   const std::vector<double>& z_U(xam_nlp->final_solution_z_U_);
   ok    &= z_L.size() == 2 && z_U.size() == 2;
   for(size_t j = 0; j < 2; ++j)
   {  ok &= 0.0 <= z_L[j] && z_L[j] <= 10. * tol;
      ok &= 0.0 <= z_U[j] && z_U[j] <= 10. * tol;
   }

   return ok;
}
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_run_xam}
*/
