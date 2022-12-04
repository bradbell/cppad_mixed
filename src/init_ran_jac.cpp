// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/is_finite_vec.hpp>
/*
{xrst_begin init_ran_jac}
{xrst_spell
   uu
}

Initialize Jacobian of Random Likelihood w.r.t. Random Effects
##############################################################

Syntax
******
*mixed_object* . ``init_ran_jac`` ( *fixed_vec* , *random_vec* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variable
:ref:`init_ran_like@init_ran_like_done_` is true.

init_ran_jac_done\_
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

   ``const CppAD::vector<a1_double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta`.

random_vec
**********
This argument has prototype

   ``const CppAD::vector<a1_double>&`` *random_vec*

It specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
{xrst_toc_hidden
   example/private/ran_jac_fun.cpp
}
ran_jac_a1fun\_
***************
The input value of the member variable

   ``CppAD::ADFun<a1_double> ran_jac_a1fun_``

does not matter.
upon return zero order forward mode for this function computes
the Jacobian of the random likelihood  with respect to the random effects;
i.e.,

.. math::

   f_u ( \theta , u )

see
:ref:`ran_jac_fun.cpp-name` .

ran_jac2hes_rc\_
****************
The input value of the member variable

   ``CppAD::mixed::sparse_rc ran_jac2hes_rc_``

does not matter.
Upon return it contains the sparsity pattern for
:math:`f_{uu} ( \theta , u )` (not just lower triangle).
The row indices are in the vector :math:`u`; i.e., just the random effects.
The column indices are in the vector  :math:`( \theta , u )`; i.e.,
both fixed and random effects.

{xrst_comment
   %example/private/ran_jac_fun_.cpp
}

{xrst_end init_ran_jac}
----------------------------------------------------------------------------
*/
void cppad_mixed::init_ran_jac(
   const d_vector& fixed_vec     ,
   const d_vector& random_vec    )
{  assert( ! init_ran_jac_done_ );
   assert( init_ran_like_done_ );
   //
   size_t m      = 1;
   size_t n_both = n_fixed_ + n_random_;
   assert( fixed_vec.size() == n_fixed_ );
   assert( random_vec.size() == n_random_ );
   assert( ran_like_fun_.Domain() == n_both );
   assert( ran_like_fun_.Range()  == m );
   // ------------------------------------------------------------------------
   // ran_jac_a1fun_
   // ------------------------------------------------------------------------
   //
   // a1_both = (fixed_vec, random_vec)
   a1_vector a1_both(n_both);
   pack(fixed_vec, random_vec, a1_both);
   //
   // Record gradient of random likelihood
   size_t abort_op_index = 0;
   bool record_compare   = false;
   CppAD::Independent(a1_both, abort_op_index, record_compare);
   //
   // a1_w = [ 1.0 ]
   a1_vector a1_w(m);
   a1_w[0]  = 1.0;
   //
   // a1_jac
   ran_like_a1fun_.Forward(0, a1_both);
   a1_vector a1_jac_both = ran_like_a1fun_.Reverse(1, a1_w);
   a1_vector a1_jac(n_random_);
   for(size_t j = 0; j < n_random_; ++j)
      a1_jac[j] = a1_jac_both[n_fixed_ + j];
   //
   if( ! CppAD::mixed::is_finite_vec( a1_jac ) )
   {  std::string error_message =
      "init_ran_like: Jacobian of ran_likelihood not finite "
      " at starting variable values";
      fatal_error(error_message);
   }
   //
   // ran_jac_a1fun_
   CppAD::ADFun<double> ran_jac_fun(a1_both, a1_jac);
# if CPPAD_MIXED_OPTIMIZE_AD_FUNCTION
   ran_jac_fun.optimize()
# endif
   ran_jac_a1fun_ = ran_jac_fun.base2ad();
   // ------------------------------------------------------------------------
   // ran_jac2hes_rc_
   // ------------------------------------------------------------------------
   //
   // ran_jac_a1fun_ computes f_u ( theta , u ), but we are only interested
   // in the sparsity patter for f_uu ( theta, u ).
   CppAD::vector<bool> select_domain(n_both), select_range(n_random_);
   for(size_t i = 0; i < n_fixed_; i++)
      select_domain[i] = false;             // false for fixed effects
   for(size_t i = 0; i < n_random_; i++)
   {  select_domain[n_fixed_ + i] = true;   // true for random effects
      select_range[i]             = true;
   }
   sparse_rc hes_pattern;
   bool transpose = false;
   ran_jac_a1fun_.subgraph_sparsity(
      select_domain, select_range, transpose, ran_jac2hes_rc_
   );
   assert(ran_jac2hes_rc_.nr() == n_random_);
   assert(ran_jac2hes_rc_.nc() == n_both);
   //
   // ------------------------------------------------------------------------
   init_ran_jac_done_ = true;
}
