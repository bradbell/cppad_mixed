// $Id:$
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
{xrst_begin init_laplace_obj}

Initialize Second Order of Approximate Objective and It's Hessian
#################################################################

Syntax
******
*mixed_object* . ``init_laplace_obj`` ( *fixed_vec* , *random_opt* )

Prototype
*********
{xrst_literal
   // BEGIN_PROTOTYPE
   // END_PROTOTYPE
}

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variable
:ref:`init_ran_like@init_ran_like_done_`
is true.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fixed_vec
*********
This specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta` at which the initialization is done.

random_opt
**********
This specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>` optimization
at which the initialization is done.
It should be the optimal value given the fixed effects
so that the Hessian w.r.t the random effects is more likely to be
positive definite.

init_laplace_obj_done\_
***********************
The input value of this member variable must be false.
Upon return it is true.

laplace_obj_fun\_, init_laplace_obj_fun_done\_
**********************************************
These have the same specification as in :ref:`init_laplace_obj_fun-name` .
This initialization is done at the *fixed_vec* value for the
fixed effects and the corresponding optimal value for the random effects.
(This Hessian of the random likelihood w.r.t. the fixed effects is
more likelihood to be positive definite at the optimal random effects.)

laplace_obj_hes\_, init_laplace_obj_hes_done\_
**********************************************
These have the same specification as in :ref:`init_laplace_obj_hes-name` .
This initialization is done at the *fixed_vec* value for the
fixed effects and the corresponding optimal value for the random effects.

{xrst_end init_laplace_obj}
*/
# include <cppad/mixed/cppad_mixed.hpp>

// BEGIN_PROTOTYPE
void cppad_mixed::init_laplace_obj(
   const d_vector&     fixed_vec        ,
   const d_vector&     random_opt       )
// END_PROTOTYPE
{  // This initialization is not called by initialize
   // so it has its own tracing.
   if( trace_init_ )
      std::cout << "Begin cppad_mixed::init_laplace_obj\n";
   //
   assert( ! init_laplace_obj_fun_done_ );
   assert( ! init_laplace_obj_hes_done_ );
   assert( ! init_laplace_obj_done_ );

   // laplace_obj_fun_
   init_laplace_obj_fun(fixed_vec, random_opt);
   assert( init_laplace_obj_fun_done_ );
   if( trace_init_ )
      std::cout << "init_laplace_obj_fun_done_\n";

   // laplace_obj_hes_
   init_laplace_obj_hes(fixed_vec, random_opt);
   assert( init_laplace_obj_hes_done_ );
   if( trace_init_ )
      std::cout << "init_laplace_obj_hes_done_\n";

   // init_laplace_obj_done_
   init_laplace_obj_done_ = true;
   if( trace_init_ )
      std::cout << "End cppad_mixed::init_laplace_obj\n";
}
