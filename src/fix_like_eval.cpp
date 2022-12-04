// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>

/*
{xrst_begin fix_like_eval}
{xrst_spell
   fabs
}

Evaluate Fixed Likelihood
#########################

Syntax
******
*vec* = *mixed_object* . ``fix_like_eval`` ( *fixed_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fix_likelihood
**************
In the special case where :ref:`fix_likelihood-name` returns the empty vector,
*vec* is also the empty vector.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta` at which :math:`g( \theta )` is evaluated.

vec
***
The return value has prototype

   ``CppAD::vector<double>`` *vec*

and is a
:ref:`problem@Negative Log-Density Vector`
corresponding to the fixed part of the negative log-likelihood
:ref:`g(theta)<theory@Fixed Likelihood, g(theta)>` .
To be specific;

:math:`g( \theta ) =`

   *vec* [0] + ``fabs`` ( *vec* [1]) + ... ``fabs`` ( *vec* [ *s* ``-1`` ])

where *s* = *vec* . ``size`` () .
{xrst_toc_hidden
   example/private/fix_like_eval.cpp
}
Example
*******
The file :ref:`fix_like_eval.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end fix_like_eval}
*/


CppAD::vector<double> cppad_mixed::fix_like_eval(const d_vector& fixed_vec)
{  assert( init_fix_like_done_ );
   if( fix_like_fun_.size_var() == 0 )
   {  // empty vector case
      return CppAD::vector<double>(0);
   }
   assert( fix_like_fun_.Domain() == n_fixed_ );
   //
   d_vector ret = fix_like_fun_.Forward(0, fixed_vec);
   if( CppAD::hasnan( ret ) ) throw CppAD::mixed::exception(
      "fix_like_eval", "resut has a nan"
   );
   return ret;
}
