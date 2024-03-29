// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin fix_like_eval.cpp dev}

fix_like_eval: Example and Test
###############################

Private
*******
This example is not part of the
:ref:`cppad_mixed public API<base_class-name>` .

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end fix_like_eval.cpp}
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>

namespace {
   using CppAD::vector;
   using CppAD::log;
   using CppAD::AD;
   using CppAD::mixed::d_sparse_rcv;
   //
   using CppAD::mixed::a1_double;

   class mixed_derived : public cppad_mixed {
   private:
      size_t                n_fixed_;
      const vector<double>& y_;
   public:
      // constructor
      mixed_derived(
         size_t                               n_fixed       ,
         size_t                               n_random      ,
         bool                                 quasi_fixed   ,
         bool                                 bool_sparsity ,
         const CppAD::mixed::d_sparse_rcv&    A_rcv         ,
         const vector<double>&                y             ) :
         cppad_mixed(
            n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
         ),
         n_fixed_(n_fixed) ,
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

         // compute this factor once
         scalar sqrt_2pi = scalar(
             CppAD::sqrt( 8.0 * CppAD::atan(1.0)
         ));

         // initialize summation
         vec[0] = scalar(0.0);

         // for each data and random effect
         for(size_t i = 0; i < y_.size(); i++)
         {  scalar mu     = theta[0] + u[i];
            scalar sigma  = theta[1];
            scalar res    = (y_[i] - mu) / sigma;

            // This is a Gaussian term, so entire density is smooth
            vec[0]  += log(sqrt_2pi * sigma) + res * res / scalar(2.0);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
      // implementation of fix_likelihood
      template <typename Vector>
      Vector template_fix_likelihood(
         const Vector& fixed_vec  )
      {  typedef typename Vector::value_type scalar;

         Vector vec(1);

         // initialize part of log-density that is smooth
         vec[0] = scalar(0.0);

         // compute these factors once
         scalar mu     = scalar(1.0);
         scalar sqrt_2 = CppAD::sqrt( scalar(2.0) );

         for(size_t j = 0; j < n_fixed_; j++)
         {
            // This is a Laplace term
            vec[0] += CppAD::log( sqrt_2 );

            // part of the density that needs absolute value
            vec.push_back(sqrt_2 * (fixed_vec[j] - mu) );
         }
         return vec;
      }
      // a1_vector version of fix_likelihood
      virtual a1_vector fix_likelihood(const a1_vector& fixed_vec)
      {  return template_fix_likelihood( fixed_vec ); }
   public:
      //
   };
}

bool fix_like_eval_xam(void)
{
   bool   ok = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();
   double sqrt_2 = CppAD::sqrt(2.0);

   size_t n_data   = 10;
   size_t n_fixed  = 2;
   size_t n_random = n_data;
   vector<double> data(n_data), fixed_vec(n_fixed), random_vec(n_random);

   fixed_vec[0] = 2.0;
   fixed_vec[1] = 0.5;
   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      random_vec[i] = double(i) / double(n_data);
   }

   // object that is derived from cppad_mixed
   bool quasi_fixed   = false;
   bool bool_sparsity = false;
   CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, data
   );
   mixed_object.initialize(fixed_vec, random_vec);

   // compute fixed negative log-density vector
   CppAD::vector<double> vec = mixed_object.fix_like_eval(fixed_vec);

   // check smooth part
   double check = CppAD::log(2.0);
   ok &= fabs( vec[0] / check - 1.0 ) <= eps;

   // check number of absolute values
   ok &= vec.size() == n_fixed + 1;

   // check argument to absolute value
   for(size_t j = 0; j < n_fixed; j++)
   {  // note that the true value is not equal to 1.0 so can deivide by check
      check = sqrt_2 * ( fixed_vec[j] - 1.0 );
      ok &= fabs( vec[1 + j] / check - 1.0 ) <= eps;
   }

   return ok;
}
// END C++
