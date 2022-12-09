// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/is_finite_vec.hpp>

/*
{xrst_begin init_hes_cross dev}

Cross Terms of Sparse Hessian w.r.t Fixed and Random Effects
############################################################

Syntax
******

| *mixed_object* . ``init_hes_cross`` (
| |tab| *fixed_vec* , *random_vec*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variable
:ref:`init_ran_like@init_ran_like_done_` is true.

init_hes_cross_done\_
*********************
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

hes_cross\_
***********
The input value of the member variable

   ``CppAD::mixed::sparse_hes_rcv hes_cross_``

does not matter.
Upon return it contains the
:ref:`sparse_hes_info-name`
for the Hessian

.. math::

   f_{u,\theta} ( \theta , u )

see :ref:`f(theta, u)<theory@Random Likelihood, f(theta, u)>`

ran_like_fun\_, ran_like_a1fun\_
================================
Either ``ran_like_fun_`` or ``ran_like_a1fun_``
can be used for the ADFun object in the
:ref:`sparse Hessian Call<sparse_hes_info@Sparse Hessian Call@f>` .

Order
*****
The results are in column major order; i.e.,

| |tab| ``hes_cross_.subset.col`` ()[ *k* ] <= ``hes_cross_.subset.col`` ()[ *k* +1]
| |tab| ``if`` ( ``hes_cross_.subset.col`` ()[ *k* ] == ``hes_cross_.subset.col`` ()[ *k* +1] )
| |tab| |tab| ``hes_cross_.subset.row`` ()[ *k* ] < ``hes_cross_.subset.row`` ()[ *k* +1]

{xrst_toc_hidden
   example/private/hes_cross.cpp
}
Example
*******
The file :ref:`hes_cross.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end init_hes_cross}
*/

void cppad_mixed::init_hes_cross(
   const d_vector& fixed_vec         ,
   const d_vector& random_vec        )
{  assert( ! init_hes_cross_done_ );
   assert( init_ran_like_done_ );
   //
   size_t m      = 1;
   size_t n_both = n_fixed_ + n_random_;
   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( ran_like_fun_.Range()  == m );
   assert( ran_like_fun_.Domain() == n_both );
   //
   // ------------------------------------------------------------------------
   // forward Jacobian sparsity for partials w.r.t fixed effects
   sparse_rc pattern_in(n_both, n_fixed_, n_fixed_);
   for(size_t k = 0; k < n_fixed_; k++)
      pattern_in.set(k, k, k);
   bool      transpose     = false;
   bool      dependency    = false;
   bool      internal_bool = bool_sparsity_;
   sparse_rc jac_pattern;
   ran_like_fun_.for_jac_sparsity(
      pattern_in, transpose, dependency, internal_bool, jac_pattern
   );
   // -----------------------------------------------------------------------
   // reverse sparsity for partial w.r.t (theta, u) of partial w.r.t theta
   CppAD::vector<bool> select_range(m);
   select_range[0] = true;
   sparse_rc hes_pattern;
   ran_like_fun_.rev_hes_sparsity(
      select_range, transpose, internal_bool, hes_pattern
   );
   // ------------------------------------------------------------------------
   // count number of elements in subset and extended sparsity pattern
   // do not need to include first n_fixed_ rows of sparsity pattern.
   size_t n_subset = 0;
   for(size_t k = 0; k < hes_pattern.nnz(); k++)
   {  if( hes_pattern.row()[k] >= n_fixed_ )
      {  assert( hes_pattern.col()[k] < n_fixed_ );
         n_subset += 1;
      }
   }
   // ------------------------------------------------------------------------
   // subset of sparstiy pattern that we are calculating
   sparse_rc subset_pattern(n_both, n_both, n_subset);
   // column major ordering
   s_vector col_major = hes_pattern.col_major();
   size_t k_subset = 0;
   for(size_t k = 0; k < hes_pattern.nnz(); k++)
   {  size_t ell = col_major[k];
      size_t r   = hes_pattern.row()[ell];
      if( r >= n_fixed_ )
      {  size_t c   = hes_pattern.col()[ell];
         assert( c < n_fixed_ );
         subset_pattern.set(k_subset++, r, c);
      }
   }
   assert( k_subset == n_subset );
   // ------------------------------------------------------------------------
   // create a d_vector containing (theta, u)
   d_vector both(n_both);
   pack(fixed_vec, random_vec, both);
   // ------------------------------------------------------------------------
   // set hes_cross_
   d_vector w(m);
   w[0] = 1.0;
   std::string coloring = "cppad.symmetric";
   hes_cross_.subset = d_sparse_rcv( subset_pattern );
   ran_like_fun_.sparse_hes(
      both,
      w,
      hes_cross_.subset,
      subset_pattern,
      coloring,
      hes_cross_.work
   );
   if( ! CppAD::mixed::is_finite_vec( hes_cross_.subset.val() ) )
   {  std::string error_message =
      "init_ran_like: Hessian of ran_likelihood cross terms not finite "
      " at starting variable values";
      fatal_error(error_message);
   }
   //
   init_hes_cross_done_ = true;
}
