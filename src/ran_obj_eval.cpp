// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>

/*
{xrst_begin ran_obj_eval}

Evaluate Laplace Approximation and Laplace Objective
####################################################

Syntax
******
*h* = *mixed_object* . ``ran_obj_eval`` ( *fixed_vec* , *random_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Purpose
*******
This routine evaluates the Laplace approximation
:ref:`h(theta, u)<theory@Objective@Laplace Approximation, h(theta, u)>` .
Note that if the random effects are optimal,
then the Laplace approximation is equal to the Laplace objective.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

ldlt_ran_hes\_
**************
It is assumed that the member variable

   ``CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_``

was updated using :ref:`update_factor-name` for the specified values of the
fixed and random effects.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta` at which :math:`h( \theta , u)` is evaluated.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

It specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u` at which :math:`h( \theta , u)` is evaluated.
Note that the Laplace approximation is equal to the Laplace objective when
:math:`u` is the
:ref:`optimal random effects<theory@Optimal Random Effects, u^(theta)>`
:math:`\hat{u} ( \theta )`.

h
*
The return value has prototype

   ``double`` *h*

and is the value of the Laplace approximation.
{xrst_toc_hidden
   example/private/ran_obj_eval.cpp
}
Example
*******
The file :ref:`ran_obj_eval.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end ran_obj_eval}
*/

// ----------------------------------------------------------------------------
double cppad_mixed::ran_obj_eval(
   const d_vector& fixed_vec  ,
   const d_vector& random_vec )
{  assert( init_ran_like_done_ );
   assert( init_ran_hes_done_ );

   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );

   // pack fixed and random effects into one vector
   d_vector both(n_fixed_ + n_random_);
   pack(fixed_vec, random_vec, both);

   // compute the logdet( f_{u,u}(theta, u )
   size_t negative;
   double logdet = ldlt_ran_hes_.logdet(negative);
   if( negative != 0 )
   {
      throw CppAD::mixed::exception(
         "ran_obj_eval",
         "The Hessian w.r.t. random effects is not positive definite."
      );
   }

   // constant term
   double pi   = CppAD::atan(1.0) * 4.0;
   double constant_term = CppAD::log(2.0 * pi) * double(n_random_) / 2.0;

   // f(theta , u)
   d_vector vec = ran_like_fun_.Forward(0, both);
   assert( vec.size() == 1);
   if( CppAD::hasnan( vec ) ) throw CppAD::mixed::exception(
      "ran_obj_eval", "result has nan"
   );

   // h(theta, u)
   double h = logdet / 2.0 + vec[0] - constant_term;

   return h;
}
