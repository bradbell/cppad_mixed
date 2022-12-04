// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin init_laplace_obj_fun}
{xrst_spell
   dyn
   ind
   nr
}

Second Order Representation of Laplace Objective and Constraints
################################################################

Syntax
******
*mixed_object* . ``init_laplace_obj_fun`` ( *fixed_vec* , *random_opt* )

Private
*******
This ``cppad_mixed`` is a :ref:`private_base_class-name` member function.

Assumptions
***********
The member variable
:ref:`init_ran_like@init_ran_like_done_`
is true.

init_laplace_obj_fun_done\_
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

It specifies the value of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u` at which the initialization is done.
It should be the optimal value given the fixed effects
so that the Hessian w.r.t the random effects is more likely to be
positive definite.

laplace_obj_fun\_
*****************
The input value of the member variable

   ``CppAD::ADFun<double> laplace_obj_fun_``

does not matter.
Upon return it contains a second order accurate recording of the
approximate Laplace objective; see
:ref:`H(beta, theta, u)<theory@Approximate Laplace Objective, H(beta, theta, u)>`
:math:`H( \beta , \theta , u )`,
followed by the approximate random constraint function; see
:ref:`B(beta, theta, u)<theory@Approximate Random Constraint Function, B(beta, theta, u)>`
:math:`B( \beta , \theta , u )`.

beta
====
The vector *beta* corresponds to the independent variables; e.g.,
``laplace_obj_fun_.Domain() == n_fixed_`` .

theta, u
========
The vector ( *theta* , *u* ) corresponds to the dynamic parameters; e.g.,
``laplace_obj_fun_.size_dyn_ind() == n_fixed_ + n_random_`` .

Range Space
===========
The function result corresponds to ( *H* , *B* ) ; e.g.,
``laplace_obj_fun_.Range() == 1 + A_rcv.nr()`` .

{xrst_end init_laplace_obj_fun}
*/
# include <Eigen/Sparse>
# include <cppad/mixed/ran_like_hes.hpp>
# include <cppad/mixed/order2random.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/is_finite_vec.hpp>

// ----------------------------------------------------------------------------
void cppad_mixed::init_laplace_obj_fun(
   const d_vector& fixed_vec  ,
   const d_vector& random_opt )
{  assert( ! init_laplace_obj_fun_done_ );
   assert( init_ran_like_done_ );
   //
   assert( A_rcv_.nnz() == A_rcv_.row().size() );
   assert( A_rcv_.nnz() == A_rcv_.col().size() );
   assert( A_rcv_.nnz() == A_rcv_.val().size() );
   //
   //  beta
   a1_vector beta(n_fixed_);
   for(size_t j = 0; j < n_fixed_; ++j)
      beta[j] = fixed_vec[j];
   //
   // theta_u
   a1_vector theta_u( n_fixed_ + n_random_ );
   pack(fixed_vec, random_opt, theta_u);
   //
   // start recording a1_double operations
   // beta:    independent variables
   // theta_u: dynamic parameters
   size_t abort_op_index = 0;
   bool record_compare   = false;
   CppAD::Independent(beta, abort_op_index, record_compare, theta_u);
   //
   // theta, u
   a1_vector theta(n_fixed_), u(n_random_);
   for(size_t j = 0; j < n_fixed_; ++j)
      theta[j] = theta_u[j];
   for(size_t j = 0; j < n_random_; ++j)
      u[j] = theta_u[j + n_fixed_];
   //
   // ran_hes_uu_rcv = f_{uu} (theta , u).
   a1_sparse_rcv ran_hes_uu_rcv;
   ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
      n_fixed_,
      n_random_,
      ran_jac_a1fun_,
      a1_ldlt_ran_hes_.pattern(),
      theta_u
   );
   //
   // a1_ldlt_ran_hes_.update
   bool ok = a1_ldlt_ran_hes_.update( ran_hes_uu_rcv );
   if( ! ok )
   {  CppAD::mixed::exception e(
         "init_laplace_obj_fun", "Hessian w.r.t. random effects is singular"
      );
      throw(e);
   }
   //
   // W(beta, theta, u)
   a1_vector W = CppAD::mixed::order2random(
      n_fixed_,
      n_random_,
      ran_jac_a1fun_,
      a1_ldlt_ran_hes_,
      beta,
      theta_u
   );
   // -----------------------------------------------------------------------
   // Evaluate f_{uu} (beta , W).
   //
   // beta_W
   a1_vector beta_W(n_fixed_ + n_random_);
   pack(beta, W, beta_W);
   //
   // ran_hes_uu_rcv
   ran_hes_uu_rcv = CppAD::mixed::ran_like_hes(
      n_fixed_,
      n_random_,
      ran_jac_a1fun_,
      a1_ldlt_ran_hes_.pattern(),
      beta_W
   );
   //
   // ldlt_obj.update
   a1_ldlt_ran_hes_.update( ran_hes_uu_rcv );
   //
   // logdet [ f_{uu} ( beta , W ) ]
   size_t negative;
   a1_double logdet = a1_ldlt_ran_hes_.logdet(negative);
   if( negative != 0 )
   {  std::string error_message =
      "init_laplace_fun: Hessian of ran_likelihood w.r.t. random effects"
      "\nnot positive at starting fixed and optimal random effects.";
      fatal_error(error_message);
   }
   //
   // Evaluate the random likelihood using (beta, W)
   a1_vector f(1);
   f  = ran_like_a1fun_.Forward(0, beta_W);
   if( ! CppAD::mixed::is_finite_vec(f) ) throw CppAD::mixed::exception(
      "init_laplace_obj", "result is not finite"
   );
   //
   // constant term
   double pi   = CppAD::atan(1.0) * 4.0;
   double constant_term = CppAD::log(2.0 * pi) * double(n_random_) / 2.0;
   //
   // now the random part of the Laplace objective
   a1_vector HB(1 + A_rcv_.nr());
   HB[0] = logdet / 2.0 + f[0] - constant_term;
   //
   if( A_rcv_.nr() > 0 )
   {
      // multiply the matrix A times the vector W and put result in
      // HB (with an index offset of 1)
      for(size_t i = 0; i < A_rcv_.nr(); i++)
         HB[1 + i] = a1_double(0.0);

      size_t K = A_rcv_.row().size();
      for(size_t k = 0; k < K; k++)
      {  size_t i = A_rcv_.row()[k];
         size_t j = A_rcv_.col()[k];
         double v = A_rcv_.val()[k];
         HB[1 + i]  += a1_double(v) * W[j];
      }
   }
   //
   laplace_obj_fun_.Dependent(beta, HB);
   laplace_obj_fun_.check_for_nan(true);
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
   std::string options =
      "no_conditional_skip no_compare_op no_print_for_op";
   laplace_obj_fun_.optimize(options);
# endif
   //
   init_laplace_obj_fun_done_ = true;
   return;
}
