// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/is_finite_vec.hpp>

/*
{xrst_begin init_ran_like}

Initialize Random Likelihood
############################

Syntax
******
*mixed_object* . ``init_ran_like`` ( *fixed_vec* , *random_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

init_ran_like_done\_
********************
The input value of this member variable must be false.
Upon return it is true.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta` at which the initialization is done.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

It specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u` at which the initialization is done.

ran_like_fun\_
**************
The input value of the member variable

   ``CppAD::ADFun<double> ran_like_fun_``

does not matter.
Upon return it contains a recording of the function
:ref:`ran_likelihood-name` .

ran_like_a1fun\_
****************
The input value of the member variable

   ``CppAD::ADFun<double> ran_like_a1fun_``

does not matter.
Upon return it contains a recording of the function
:ref:`ran_likelihood-name` .

{xrst_end init_ran_like}
*/

void cppad_mixed::init_ran_like(
   const d_vector& fixed_vec  ,
   const d_vector& random_vec )
{  assert( ! init_ran_like_done_ );
   //
   using CppAD::AD;
   using CppAD::ADFun;
   using CppAD::vector;
   using CppAD::Independent;
   using CppAD::mixed::is_finite_vec;
   //
   // ------------------------------------------------------------------
   // record ran_like_fun_
   // ------------------------------------------------------------------
   // combine into one vector
   a1_vector a1_both( n_fixed_ + n_random_ );
   pack(fixed_vec, random_vec, a1_both);

   // start recording a1_double operations
   size_t abort_op_index = 0;
   bool record_compare   = false;
   Independent(a1_both, abort_op_index, record_compare);

   // extract the fixed and random effects
   a1_vector a1_theta(n_fixed_), a1_u(n_random_);
   unpack(a1_theta, a1_u, a1_both);

   // compute ran_likelihood using a1_double operations
   a1_vector a1_vec = ran_likelihood(a1_theta, a1_u);
   if( a1_vec.size() == 0 )
   {  std::string error_message =
         "init_ran_like: n_random > 0 and ran_likelihood has size 0";
      fatal_error(error_message);
   }
   if( a1_vec.size() != 1 )
   {  std::string error_message =
      "init_ran_like: ran_likelihood does not have size zero or one.";
      fatal_error(error_message);
   }
   if( ! is_finite_vec( a1_vec ) )
   {  std::string error_message =
      "init_ran_like: ran_likelihood not finite at starting variable values";
      fatal_error(error_message);
   }

   // save the recording
   ran_like_fun_.Dependent(a1_both, a1_vec);
   ran_like_fun_.check_for_nan(true);

   // optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
   std::string options =
      "no_conditional_skip no_compare_op no_print_for_op";
   ran_like_fun_.optimize(options);
# endif
   // ------------------------------------------------------------------
   // set ran_like_a1fun_
   // ------------------------------------------------------------------
   ran_like_a1fun_ = ran_like_fun_.base2ad();
   //
   // ------------------------------------------------------------------
   init_ran_like_done_ = true;
   return;
}
