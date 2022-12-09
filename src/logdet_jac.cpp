// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
{xrst_begin logdet_jac dev}
{xrst_spell
   determinant
   logdet
}

Jacobian of Log Determinant of Hessian w.r.t. Random Effects
############################################################

Syntax
******

| *mixed_object* . ``logdet_jac`` (
| |tab| *fixed_vec* , *random_vec* , *logdet_fix* , *logdet_ran*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Purpose
*******
This routine computes the Jacobian of the log determinant
of the Hessian of the random likelihood
:ref:`f(theta, u)<theory@Random Likelihood, f(theta, u)>`
with respect to the random effects vector :math:`u`.
To be specific, it computes both

.. math::

   \partial_\theta \log \det [ f_{u,u} ( \theta, u ) ]
   \; \R{and} \;
   \partial_u \log \det [ f_{u,u} ( \theta, u ) ]

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

ldlt_ran_hes\_
**************
It is assumed that the member variable

   ``CPPAD_MIXED_LDLT ldlt_ran_hes_``

was updated using :ref:`update_factor-name` for the specified values of the
fixed and random effects.

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

and is the value of fixed effects :math:`\theta`.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_vec*

and is the value of fixed effects :math:`u`.

logdet_fix
**********
This argument has prototype

   ``CppAD::vector<double>&`` *logdet_fix*

Its input size must be equal to ``n_fixed_`` .
Upon return, it contains the value of the derivative w.r.t
the fixed effects.
``ran_hes_.col`` [ *k* ] ``- n_fixed_`` of the Hessian.

logdet_ran
**********
This argument has prototype

   ``CppAD::vector<double>&`` *logdet_ran*

Its input size must be equal to ``n_random_`` .
Upon return, it contains the value of the derivative w.r.t
the random effects.
{xrst_toc_hidden
   example/private/logdet_jac.cpp
}
Example
*******
The file :ref:`logdet_jac.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end logdet_jac}
*/
// ----------------------------------------------------------------------------
void cppad_mixed::logdet_jac(
   const d_vector& fixed_vec  ,
   const d_vector& random_vec ,
   d_vector&       logdet_fix ,
   d_vector&       logdet_ran )
{  assert( init_ran_hes_done_ );
   //
   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( logdet_fix.size() == n_fixed_ );
   assert( logdet_ran.size() == n_random_ );
   //
   // compute the inverse where Hessian is possibly non-zero
   CppAD::mixed::sparse_mat_info weight_info;
   size_t K = ran_hes_uu_rcv_.nnz();
   weight_info.row.resize(K);
   weight_info.col.resize(K);
   weight_info.val.resize(K);
   for(size_t k = 0; k < K; k++)
   {  size_t r = ran_hes_uu_rcv_.row()[k];
      size_t c = ran_hes_uu_rcv_.col()[k];
      assert(r < n_random_);
      assert(c < n_random_);
      weight_info.row[k] = r;
      weight_info.col[k] = c;
   }
   ldlt_ran_hes_.inv(
         weight_info.row,
         weight_info.col,
         weight_info.val
   );
   // must weight the off diagonal elements twice
   // (to account for upper diagonal entry which is not present)
   for(size_t k = 0; k < K; k++)
   {  if( weight_info.row[k] != weight_info.col[k] )
         weight_info.val[k] *= 2.0;
   }
   d_vector dw(n_fixed_ + n_random_);
   dw = ran_hes_fun_.Reverse(1, weight_info.val);
   if( CppAD::hasnan( dw ) ) throw CppAD::mixed::exception(
      "logdet_jac", "result has a nan"
   );
   //
   // split out fixed and random parts of the derivative
   unpack(logdet_fix, logdet_ran, dw);
   //
   return;
}
