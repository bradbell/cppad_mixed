// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
{xrst_begin init_laplace_obj_hes dev}

Initialize Hessian of Approximate Laplace Objective
###################################################

Syntax
******

| *mixed_object* . ``init_laplace_obj_hes`` (
| |tab| *fixed_vec* , *random_opt*
| )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variable
:ref:`init_laplace_obj_fun@init_laplace_obj_fun_done_`
is true.

init_laplace_obj_hes_done\_
***************************
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

random_opt
**********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_opt*

It specifies the initial value for the
:ref:`random effects<problem@Notation@Random Effects, u>` optimization.
It should be the optimal value given the fixed effects
so that the Hessian w.r.t the random effects is more likely to be
positive definite.

laplace_obj_hes\_
*****************
The input value of the member variable

   ``CppAD::mixed::sparse_hes_info laplace_obj_hes_``

does not matter.
Upon return it contains the
:ref:`sparse_hes_info-name`
for the lower triangle of the Hessian

.. math::

   r_{\theta,\theta} ( \theta )
   =
   H_{\beta,\beta} [ \beta, \theta , \hat{u} ( \theta) ]

see
:ref:`H(beta, theta, u)<theory@Approximate Laplace Objective, H(beta, theta, u)>` .
Note that the matrix is symmetric and hence can be recovered from
its lower triangle.

laplace_obj_fun\_
=================
This ADFun object should be used in the
:ref:`sparse Hessian Call<sparse_hes_info@Sparse Hessian Call@f>` .

{xrst_end init_laplace_obj_hes}
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_laplace_obj_hes(
   const d_vector& fixed_vec     ,
   const d_vector& random_opt    )
{  assert( ! init_laplace_obj_hes_done_ );
   assert( init_laplace_obj_fun_done_ );
   //
   // beta
   d_vector beta = fixed_vec;
   //
   // theta_u
   d_vector theta_u(n_fixed_ + n_random_);
   pack(fixed_vec, random_opt, theta_u);
   //
   // Compute Jacobian sparsity for partials w.r.t. beta.
   // Note that theta and u are dynamic parameters in laplace_obj_fun_.
   sparse_rc pattern_in(n_fixed_, n_fixed_, n_fixed_);
   for(size_t k = 0; k < n_fixed_; k++)
      pattern_in.set(k, k, k);
   bool      transpose     = false;
   bool      dependency    = false;
   bool      internal_bool = bool_sparsity_;
   sparse_rc jac_pattern;
   laplace_obj_fun_.for_jac_sparsity(
      pattern_in, transpose, dependency, internal_bool, jac_pattern
   );

   // compute sparsity pattern corresponding to
   // H_{beta,beta} (beta, theta, u)
   // Note that theta and u are dynamic parameters in laplace_obj_fun_.
   size_t m = laplace_obj_fun_.Range();
   CppAD::vector<bool> select_range(m);
   for(size_t i = 0; i < m; i++)
      select_range[i] = false;
   select_range[0] = true;
   sparse_rc hes_pattern;
   laplace_obj_fun_.rev_hes_sparsity(
      select_range, transpose, internal_bool, hes_pattern
   );
   assert( hes_pattern.nr() == n_fixed_ );
   assert( hes_pattern.nc() == n_fixed_ );
   //
   // beta_beta_pattern
   // sparsity pattern for lower traingle of H_{beta,beta} (beta, theta, u)
   size_t nnz = 0;
   for(size_t k = 0; k < hes_pattern.nnz(); k++)
   {  size_t r = hes_pattern.row()[k];
      size_t c = hes_pattern.col()[k];
      if( r >= c )
         ++nnz;
   }
   sparse_rc beta_beta_pattern(n_fixed_, n_fixed_, nnz);
   size_t ell = 0;
   for(size_t k = 0; k < hes_pattern.nnz(); k++)
   {  size_t r = hes_pattern.row()[k];
      size_t c = hes_pattern.col()[k];
      if( r >= c )
         beta_beta_pattern.set(ell++, r, c);
   }
   assert( nnz == ell );
   //
   // laplace_obj_hes_.subset
   // only compute the lower triangle of H_beta_beta ( beta, theta, u)
   laplace_obj_hes_.subset = d_sparse_rcv( beta_beta_pattern );
   //
   // compute work for reuse during sparse Hessian calculations
   // (value of weight vecgttor does not affectg work, but avoid wanrings)
   d_vector weight(m);
   for(size_t i = 0; i < m; i++)
      weight[i] = 1.0;
   std::string coloring  = "cppad.symmetric";
   laplace_obj_fun_.sparse_hes(
      beta                   ,
      weight                 ,
      laplace_obj_hes_.subset ,
      beta_beta_pattern      ,
      coloring               ,
      laplace_obj_hes_.work
   );
   //
   init_laplace_obj_hes_done_ = true;
}
