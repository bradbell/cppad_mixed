// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_random.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_random_ctor}
{xrst_spell
   nlp
   nnz
}

Ipopt Random Optimization Callback Constructor
##############################################

Syntax
******

| ``CppAD::mixed::ipopt_random`` *ipopt_object* (
| |tab| *fixed_vec* ,
| |tab| *random_lower* ,
| |tab| *random_upper* ,
| |tab| *random_in* ,
| |tab| *mixed_object*
| )

fixed_vec
*********
specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta`.
It is stored as a reference so it must exist for as long as
*ipopt_object* exists.

random_lower
************
this vector has size
:ref:`derived_ctor@n_random`
and specifies the lower limits for the optimization of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
It is stored as a reference so it must exist for as long as
*ipopt_object* exists.
The value minus infinity can be used to specify no lower limit.

random_upper
************
this vector has size *n_random*
and specifies the upper limits for the optimization of the random effects.
It is stored as a reference so it must exist for as long as
*ipopt_object* exists.
The value plus infinity can be used to specify no upper limit.

random_in
*********
this vector has size *n_random* and specifies
the initial value used for the optimization of the random effects.
It is stored as a reference so it must exist for as long as
*ipopt_object* exists.
The value plus infinity can be used to specify no upper limit.

mixed_object
************
The argument *mixed_object* is an object of a class that is
derived from the ``cppad_mixed`` base class.

Member Variables
****************

nlp_lower_bound_inf\_
=====================
set to a finite value that is used by Ipopt for minus infinity.

nlp_upper_bound_inf\_
=====================
set to a finite value that is used by Ipopt for plus infinity.

nnz_h_lag\_
===========
set to the number of non-zero entries in the Hessian of the Lagrangian.
This is the same as for the Hessian of the objective because there
are no constraints (except for box constraints) in this problem.

{xrst_spell_off}
{xrst_code cpp} */
ipopt_random::ipopt_random(
   const d_vector&     fixed_vec          ,
   const d_vector&     random_lower       ,
   const d_vector&     random_upper       ,
   const d_vector&     random_in          ,
   cppad_mixed&        mixed_object       ) :
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_random_ctor}
*/
n_fixed_             ( fixed_vec.size()  )                ,
n_random_            ( random_lower.size() )              ,
fixed_vec_           ( fixed_vec )                        ,
random_lower_        ( random_lower )                     ,
random_upper_        ( random_upper )                     ,
random_in_           ( random_in )                        ,
nnz_h_lag_           ( mixed_object.ran_hes_uu_rcv_.nnz() )  ,
mixed_object_        ( mixed_object    )                  ,
error_random_        ( n_random_ )
{  // -----------------------------------------------------------------------
   // set nlp_lower_bound_inf_, nlp_upper_bound_inf_
   double inf           = std::numeric_limits<double>::infinity();
   nlp_lower_bound_inf_ = - 1e19;
   nlp_upper_bound_inf_ = + 1e19;
   for(size_t j = 0; j < n_random_; j++)
   {  if( random_lower[j] != - inf ) nlp_lower_bound_inf_ =
            std::min(nlp_lower_bound_inf_, 1.1 * random_lower[j] );
      //
      if( random_upper[j] != inf ) nlp_upper_bound_inf_ =
            std::max(nlp_upper_bound_inf_, 1.1 * random_upper[j] );
   }
   objective_current_ = std::numeric_limits<double>::quiet_NaN();
   return;
}
} } // END_CPPAD_MIXED_NAMESPACE
