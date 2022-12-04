// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
{xrst_begin update_factor}

Update the Factorization of Hessian w.r.t. Random Effects
#########################################################

Syntax
******
*mixed_object* . ``update_factor`` ( *fixed_vec* , *random_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Purpose
*******
This routine updates the Cholesky factorization of the Hessian
:math:`f_{u,u} ( \theta , u )^{-1}`
so that it corresponds to the current value of
:math:`\theta` and :math:`u`.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

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
If this factorization is used for computing derivatives of the
:ref:`Laplace objective<theory@Objective@Laplace Objective, r(theta)>` ,
this should be the
:ref:`optimal random effects<theory@Optimal Random Effects, u^(theta)>` .

ran_hes_fun\_
*************
The :ref:`private_base_class@ran_hes_fun_` member variable
will hold the first order Taylor coefficient corresponding
to the specified fixed and random effects; i.e.,

   ``ran_hes_fun_.Forward`` (0, *both* )

has been called where *both* is a packed version
of the fixed and random effects.

ldlt_ran_hes\_
**************
On input, the member variable

   ``CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_``

has been
:ref:`initialized<ldlt_eigen_init-name>`
using the sparsity pattern for the Hessian.
Upon return, ``ldlt_ran_hes_`` contains the updated
:ref:`factorization<ldlt_eigen_update-name>`
corresponding to the specified values for the fixed
and random effects.
{xrst_toc_hidden
   example/private/update_factor.cpp
}
Example
*******
The file :ref:`update_factor.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end update_factor}
*/
// ----------------------------------------------------------------------------
void cppad_mixed::update_factor(
   const d_vector& fixed_vec  ,
   const d_vector& random_vec )
{  assert( init_ran_hes_done_ );
   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   //
   // pack fixed and random effects into one vector
   d_vector both(n_fixed_ + n_random_);
   pack(fixed_vec, random_vec, both);
   //
   // set the value vector in the sparse matrix information
   d_vector hes_val = ran_hes_fun_.Forward(0, both);
   if( CppAD::hasnan( hes_val ) ) throw CppAD::mixed::exception(
      "update_factor", "result has nan"
   );
   //
   // ran_hes_uu_rcv_
   size_t nnz = ran_hes_uu_rcv_.nnz();
   assert( hes_val.size() == nnz );
   for(size_t k = 0; k < nnz; ++k)
      ran_hes_uu_rcv_.set(k, hes_val[k]);
   //
   // update the LDLT factor
   bool ok = ldlt_ran_hes_.update(ran_hes_uu_rcv_);
   if( ! ok )
   {  CppAD::mixed::exception e(
         "update_factor", "Hessian w.r.t. random effects is singular"
      );
      throw(e);
   }
   return;
}
