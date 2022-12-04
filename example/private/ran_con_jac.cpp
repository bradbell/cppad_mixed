// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin ran_con_jac.cpp}

ran_con_jac: Example and Test
#############################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

Model
*****

.. math::

   \B{p}( y_i | \theta , u ) \sim \B{N} ( u_i + \theta_0 , \theta_1^2 )

.. math::

   \B{p}( u_i | \theta ) \sim \B{N} ( 0 , 1 )

.. math::

   A = [ 1 , \cdots , 1 ]

It follows that the *i*-th component of the
optimal random effects minimize the function

.. math::

   0.5 ( u_i + theta_0 - y_i )^2 / \theta_1^2 + 0.5 u_i^2

Taking the derivative w.r.t :math:`u_i` and setting it equal to zero
we have

.. math::

   \hat{u}_i ( \theta ) = ( y_i - \theta_0 ) / ( 1 + \theta_1^2 )

The random constraint function is

.. math::

   \sum_i \hat{u}_i ( \theta )
   =
   \sum_i
   ( y_i - \theta_0 ) / ( 1 + \theta_1^2 )

The partial w.r.t. :math:`\theta_0` is

.. math::

   - \sum_i  1.0 / ( 1 + \theta_1^2 )

The partial w.r.t. :math:`\theta_1` is

.. math::

   - 2 \sum_i ( y_i - \theta_0 ) \theta_1 / ( 1 + \theta_1^2 )

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end ran_con_jac.cpp}
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
   using CppAD::vector;
   using CppAD::log;
   using CppAD::AD;
   using CppAD::mixed::sparse_mat_info;
   //
   using CppAD::mixed::a1_double;

   class mixed_derived : public cppad_mixed {
   private:
      const vector<double>& y_;
   public:
      // constructor
      mixed_derived(
         size_t                               n_fixed       ,
         size_t                               n_random      ,
         bool                                 quasi_fixed   ,
         bool                                 bool_sparsity ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv          ,
         const vector<double>&                y             ) :
         cppad_mixed(
            n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
         ),
         y_(y)
      {  assert( n_fixed == 2);
      }
      // implementation of ran_likelihood
      template <typename Vector>
      Vector template_ran_likelihood(
         const Vector& theta  ,
         const Vector& u      )
      {  typedef typename Vector::value_type scalar;

         Vector vec(1);

         // initialize part of log-density that is always smooth
         vec[0] = scalar(0.0);

         // pi
         scalar sqrt_2pi = scalar(
             CppAD::sqrt(8.0 * CppAD::atan(1.0)
         ));

         for(size_t i = 0; i < y_.size(); i++)
         {  scalar mu     = u[i] + theta[0];
            scalar sigma  = theta[1];
            scalar res    = (y_[i] - mu) / sigma;

            // p(y_i | u, theta)
            vec[0] += log(sqrt_2pi * sigma) + res*res / scalar(2.0);

            // p(u_i | theta)
            vec[0] += log(sqrt_2pi) + u[i] * u[i] / scalar(2.0);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
   public:
      //
   };
}

bool ran_con_jac_xam(void)
{
   bool   ok = true;

   size_t n_data   = 10;
   size_t n_fixed  = 2;
   size_t n_random = n_data;
   vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);
   vector<double> uhat(n_random);

   fixed_vec[0] = 2.0;
   fixed_vec[1] = 1.0;
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 2);
      random_vec[i] = 0.0;
   }

   // lower and upper limits for random effects
   double inf = std::numeric_limits<double>::infinity();
   vector<double> random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }

   // constraint matrix will sum all the random effects
   // nr = 1, nc = n_random, nnz = n_random
   CppAD::mixed::sparse_rc A_pattern(1, n_random, n_random);
   for(size_t j = 0; j < n_random; j++)
      A_pattern.set(j, 0, j);
   CppAD::mixed::d_sparse_rcv A_rcv(A_pattern);
   for(size_t j = 0; j < n_random; j++)
      A_rcv.set(j, 1.0);

   // mixed_object
   bool quasi_fixed   = false;
   bool bool_sparsity = true;
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
   );
   mixed_object.initialize(fixed_vec, random_vec);

   // optimize the random effects
   std::string options;
   options += "Integer print_level 0\n";
   options += "Numeric tol         1e-10\n";
   options += "String  sb          yes\n";
   options += "String  derivative_test second-order\n";
   uhat = mixed_object.optimize_random(
      options, fixed_vec, random_lower, random_upper, random_vec
   );

   // must factor f_{u,u} (theta, uhat)
   mixed_object.update_factor(fixed_vec, uhat);

   // compute sparstiy pattern for jacobian of random constraints
   CppAD::mixed::d_sparse_rcv jac_rcv;
   mixed_object.ran_con_jac(fixed_vec, uhat, jac_rcv);

   // check number of possibly non_zero elements.
   ok &= jac_rcv.nnz() == 2;
   //
   // partial w.r.t. theta_0
   ok &= jac_rcv.row()[0] == 0;
   ok &= jac_rcv.col()[0] == 0;
   //
   // partial w.r.t. theta_1
   ok &= jac_rcv.row()[1] == 0;
   ok &= jac_rcv.col()[1] == 1;

   // Now compute the sparse Jacobian values
   mixed_object.ran_con_jac(fixed_vec, uhat, jac_rcv);
   //
   double jac_0 = jac_rcv.val()[0];
   double jac_1 = jac_rcv.val()[1];

   // check results
   double check_0 = 0.0;
   double check_1 = 0.0;
   for(size_t i = 0; i < n_random; i++)
   {  double theta_0   = fixed_vec[0];
      double theta_1   = fixed_vec[1];
      double theta_1sq = theta_1 * theta_1;
      check_0 -= 1.0 / (1.0 + theta_1sq );
      check_1 -= (data[i] - theta_0) * theta_1 / (1.0 + theta_1sq );
   }

   // used 1e-10 for random_optmize tolerance
   ok &= fabs( jac_0 / check_0 - 1.0 ) < 1e-9;
   ok &= fabs( jac_1 / check_1 - 1.0 ) < 1e-9;

   return ok;
}
// END C++
