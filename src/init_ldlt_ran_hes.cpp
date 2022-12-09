// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin init_ldlt_ran_hes dev}
{xrst_spell
   uu
}

Initialize Cholesky Factor of Hessian of Random Likelihood
##########################################################

Syntax
******
*mixed_object* . ``init_ldlt_ran_hes`` ()

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variable
:ref:`init_ran_hes@init_ran_hes_done_` is true.

init_ldlt_ran_hes_done\_
************************
The input value of this member variable must be false.
Upon return it is true.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

ldlt_ran_hes\_
**************
The member variable

   ``CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_``

must not been previously initialized

   ``size_t ldlt_ran_hes_``

Upon return, the function

   ``ldlt_ran_hes_.init`` ( ``ran_hes_uu_rcv_.pat`` ())

has been called with *ran_hes_uu_rcv_.pat* ()
equal to the sparsity pattern for the
Hessian of the random likelihood with respect to the random effects;
see :ref:`ldlt_eigen_init-name` .
This sparsity information is relative to just the random effects;
i.e., the row and column indices are relative to the vector :math:`u )`.
For each row and column indices in *ran_hes_uu_rcv_.pat* () ,
the corresponding row and column index in ``ran_hes_uu_rcv_`` is
``n_fixed_`` greater.

a1_ldlt_ran_hes\_
*****************
The member variable

   ``CppAD::mixed::ldlt_eigen<a1_double> a1_ldlt_ran_hes_``

is initialized the same as ``ldlt_ran_hes_`` .

{xrst_end init_ldlt_ran_hes}
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_ldlt_ran_hes(void)
{  assert( ! init_ldlt_ran_hes_done_ );
   assert( init_ran_hes_done_ );
   //
   //
   // initialize ldlt_ran_hes_
   ldlt_ran_hes_.init(ran_hes_uu_rcv_.pat());
   //
   // initialize a1_ldlt_ran_hes_
   a1_ldlt_ran_hes_.init(ran_hes_uu_rcv_.pat());
   //
   init_ldlt_ran_hes_done_ = true;
   return;
}
