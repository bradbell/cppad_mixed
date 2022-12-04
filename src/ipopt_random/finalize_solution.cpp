// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_finalize_solution}
{xrst_spell
   cq
   infeasibility
   infeasible
   iterates
   namespace
   naninf
   solution solution
   unrecoverable
}

Get Solution Results
####################

Syntax
******

| ``finalize_solution`` (
| |tab| *status* , *n* , *x* , *z_L* , *z_U* , *m* , *g* , *lambda* , *obj_value* , *ip_data* , *ip_cq*
| )

solution\_
**********
This routine checks the solution values and sets the member variable

   ``CppAD::mixed fixed_solution solution_``

.

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

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_finalize_solution}
*/
{  assert( size_t(n) == n_random_ );
   assert( m == 0 );
   //
   assert( random_opt_.size() == 0 );
   random_opt_.resize(n_random_);
   for(size_t j = 0; j < n_random_; j++)
      random_opt_[j] = x[j];
   //
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
