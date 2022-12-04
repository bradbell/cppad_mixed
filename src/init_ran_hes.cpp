// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ran_like_hes.hpp>
# include <cppad/mixed/is_finite_vec.hpp>

/*
{xrst_begin init_ran_hes}
{xrst_spell
   uu
}

Initialize Hessian of Random Likelihood w.r.t Random Effects
############################################################

Syntax
******

| *mixed_object* . ``init_ran_hes`` (
| |tab| *fixed_vec* , *random_vec*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variables
:ref:`init_ran_like@init_ran_like_done_` and
:ref:`init_ran_jac@init_ran_jac_done_` are true.

init_ran_hes_done\_
*******************
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

ran_hes_uu_rcv\_
****************
The input value of the member variable

   ``CppAD::mixed::d_sparse_rcv ran_hes_uu_rcv_``

does not matter.
Upon return ``ran_hes_uu_rcv_.pat()`` contains the sparsity pattern
for the lower triangle of the Hessian

.. math::

   f_{u,u} ( \theta , u )

see :ref:`f(theta, u)<theory@Random Likelihood, f(theta, u)>`
The matrix is symmetric and hence can be recovered from
its lower triangle.

Random Effects Index
====================
The indices in *ran_hes_uu_rcv_* are relative to just
the random effects and hence are all less that ``n_random_`` .

Order
=====
The results are in column major order; i.e.,

| |tab| ``ran_hes_uu_rcv_.col`` ()[ *k* ] <= ``ran_hes_uu_rcv_.col`` ()[ *k* +1]
| |tab| ``if`` ( ``ran_hes_uu_rcv_.col`` ()[ *k* ] == ``ran_hes_uu_rcv_.col`` ()[ *k* +1] )
| |tab| |tab| ``ran_hes_uu_rcv_.row`` ()[ *k* ] < ``ran_hes_uu_rcv_.row`` ()[ *k* +1]

ran_hes_fun\_
*************
The input value of the member variables

   ``CppAD::ADFun<double> ran_hes_fun_``

does not matter.
Upon return its zero order forward mode computes
the lower triangle of the sparse Hessian

.. math::

   f_{u,u} ( \theta , u )

in the same order as the elements of
``ran_hes_uu_rcv_`` (and ``a1_hes_rcv_`` ).

Contents
********
{xrst_toc_list
   example/private/ran_hes_fun.cpp
}

{xrst_end init_ran_hes}
*/

void cppad_mixed::init_ran_hes(
   const d_vector& fixed_vec     ,
   const d_vector& random_vec    )
{  assert( ! init_ran_hes_done_ );
   assert( init_ran_like_done_ );
   assert( init_ran_jac_done_ );
   //
   size_t n_both = n_fixed_ + n_random_;
   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( ran_like_fun_.Domain() == n_both );
   assert( ran_like_fun_.Range()  == 1 );
   //
   // a1_both = (fixed_vec, random_vec)
   a1_vector a1_both(n_both);
   pack(fixed_vec, random_vec, a1_both);
   //
   // hes pattern relative to both fixed and random effects
   // (also count number of entries in lower traingle)
   size_t nnz   = ran_jac2hes_rc_.nnz();
   size_t n_low = 0;
   sparse_rc hes_pattern(n_both, n_both, nnz);
   for(size_t k = 0; k < nnz; ++k)
   {  size_t r = ran_jac2hes_rc_.row()[k];
      size_t c = ran_jac2hes_rc_.col()[k];
      assert( r < n_random_ );
      r = r + n_fixed_;
      hes_pattern.set(k, r, c);
      if( r >= c )
         ++n_low;
   }
   //
   // subset of sparstiy pattern that we are calculating
   // in column major order
   sparse_rc ran_hes_uu_rc(n_random_,  n_random_, n_low);
   s_vector col_major = hes_pattern.col_major();
   size_t k_low = 0;
   for(size_t k = 0; k < nnz; k++)
   {  size_t ell = col_major[k];
      size_t r   = hes_pattern.row()[ell];
      size_t c   = hes_pattern.col()[ell];
      assert( r >= n_fixed_ );
      assert( c >= n_fixed_ );
      if( r >= c )
      {  ran_hes_uu_rc.set(k_low, r - n_fixed_, c - n_fixed_);
         ++k_low;
      }
   }
   assert( k_low == n_low );
   //
   // ran_hes_uu_rcv_
   ran_hes_uu_rcv_ = d_sparse_rcv( ran_hes_uu_rc );
   //
   // Declare the independent and dependent variables for taping calculation
   // of Hessian of the random likelihood w.r.t. the random effects
   size_t abort_op_index = 0;
   bool record_compare   = false;
   CppAD::Independent(a1_both, abort_op_index, record_compare);
   //
   // a1_val
   a1_sparse_rcv a1_ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
      n_fixed_, n_random_, ran_jac_a1fun_, ran_hes_uu_rc, a1_both
   );
   const a1_vector& a1_val( a1_ran_hes_uu_rcv.val() );
   //
   if( ! CppAD::mixed::is_finite_vec( a1_val ) )
   {  std::string error_message =
      "init_ran_like: Hessian of ran_likelihood w.r.t random effects"
      " not finite at starting variable values";
      fatal_error(error_message);
   }
   //
   // ran_hes_fun_
   ran_hes_fun_.Dependent(a1_both, a1_val);
   ran_hes_fun_.check_for_nan(true);
   //
   // optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
   std::string options =
      "no_conditional_skip no_compare_op no_print_for_op";
   ran_hes_fun_.optimize(options);
# endif
   //
   init_ran_hes_done_ = true;
}
