// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ran_con_eval dev}

Evaluate the Random Constraint Function
#######################################

Syntax
******

   *mixed_object* . ``ran_con_eval`` ( *random_vec* , *Au* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

It specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u` at which *A* * *u* is evaluated.

Au
**
This argument has prototype

   ``CppAD::vector<double>`` *Au*

Its size must be equal to the number of rows in the
:ref:`random constraint matrix<problem@Notation@Random Constraint Matrix, A>` .
The input value of its elements does not matter.
Upon return, it contains the product *A* * *u* .
If the argument *random_vec* is the
:ref:`optimal random effects<theory@Optimal Random Effects, u^(theta)>`
*Au* is the value of the
:ref:`random constraint Function<problem@Notation@Random Constraint Function, A*u^(theta)>` .
{xrst_toc_hidden
   example/private/ran_con_eval.cpp
}
Example
*******
The file :ref:`ran_con_eval.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end ran_con_eval}
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::ran_con_eval(
   const d_vector& random_vec ,
   d_vector&       Au         )
{  assert( random_vec.size() == n_random_ );
   assert( Au.size() == A_rcv_.nr() );
   assert( A_rcv_.row().size() == A_rcv_.col().size() );
   assert( A_rcv_.row().size() == A_rcv_.val().size() );
   size_t K = A_rcv_.row().size();
   //
   for(size_t i = 0; i < A_rcv_.nr(); i++)
      Au[i] = 0.0;
   for(size_t k = 0; k < K; k++)
   {  size_t i = A_rcv_.row()[k];
      size_t j = A_rcv_.col()[k];
      double v = A_rcv_.val()[k];
      Au[i] += v * random_vec[j];
   }
   return;
}
