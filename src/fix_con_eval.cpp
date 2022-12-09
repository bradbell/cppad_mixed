// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>

/*
{xrst_begin fix_con_eval dev}

Evaluate Fixed Constraint Function
##################################

Syntax
******
*vec* = *mixed_object* . ``fix_con_eval`` ( *fixed_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

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
vector :math:`\theta` at which :math:`c( \theta )` is evaluated.

vec
***
The return value has prototype

   ``CppAD::vector<double>`` *vec*

and is the constraint function value
corresponding to the fixed effects; see
:ref:`fix_constraint@vec` .
{xrst_toc_hidden
   example/private/fix_con_eval.cpp
}
Example
*******
The file :ref:`fix_con_eval.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end fix_con_eval}
*/


CppAD::vector<double> cppad_mixed::fix_con_eval(const d_vector& fixed_vec)
{
   // make sure initialize has been called
   if( ! initialize_done_ )
   {  std::string error_message =
      "fix_con_eval: initialize was not called before constraint_eval";
      fatal_error(error_message);
   }
   if( fix_con_fun_.size_var() == 0 )
   {  return CppAD::vector<double>(0); // empty vector
   }
   assert( fix_con_fun_.Domain() == n_fixed_ );
   //
   d_vector ret = fix_con_fun_.Forward(0, fixed_vec);
   if( CppAD::hasnan( ret ) ) throw CppAD::mixed::exception(
      "fix_con_eval", "resut has a nan"
   );
   return ret;
}
