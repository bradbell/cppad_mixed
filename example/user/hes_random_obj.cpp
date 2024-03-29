// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin hes_random_obj.cpp}

Hessian of Random Effects Objective: Example and Test
#####################################################

{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end hes_random_obj.cpp}
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>


namespace {
   using CppAD::log;
   using CppAD::AD;
   //
   using CppAD::mixed::d_sparse_rcv;
   using CppAD::mixed::d_vector;
   using CppAD::mixed::a1_vector;
   using CppAD::mixed::s_vector;
   //
   class mixed_derived : public cppad_mixed {
   private:
      const d_vector&       y_;
   public:
      // constructor
      mixed_derived(
         size_t                 n_fixed       ,
         size_t                 n_random      ,
         bool                   quasi_fixed   ,
         bool                   bool_sparsity ,
         const d_vector&       y              ) :
         cppad_mixed(n_fixed, n_random, quasi_fixed, bool_sparsity) ,
         y_(y)
      { }
      // implementation of ran_likelihood
      a1_vector ran_likelihood(
         const a1_vector&         theta  ,
         const a1_vector&         u      ) override
      {
         a1_vector vec(1);

         // compute this factor once
         a1_double sqrt_2pi =  CppAD::sqrt( 8.0 * CppAD::atan(1.0) );

         // initialize summation
         vec[0] = 0.0;

         // for each data and random effect
         for(size_t i = 0; i < y_.size(); i++)
         {  a1_double mu     = u[i];
            a1_double sigma  = theta[i] + double(i) / double( y_.size() );
            a1_double res    = (y_[i] - mu) / sigma;

            // This is a Gaussian term, so entire density is smooth
            vec[0]  += log(sqrt_2pi * sigma) + res * res / 2.0;
         }
         return vec;
      }
   };
}

bool hes_random_obj_xam(void)
{
   bool   ok  = true;
   double eps = 100. * std::numeric_limits<double>::epsilon();
   //
   size_t n_data   = 10;
   size_t n_fixed  = n_data;
   size_t n_random = n_data;
   d_vector  data(n_data);
   d_vector  fixed_vec(n_fixed), random_vec(n_random);
   a1_vector a1_fixed(n_fixed), a1_random(n_random);

   for(size_t i = 0; i < n_data; i++)
   {  data[i]       = double(i + 1);
      //
      fixed_vec[i]  = 1.5;
      a1_fixed[i]   = fixed_vec[i];
      //
      random_vec[i] = 0.0;
      a1_random[i]  = random_vec[i];
   }

   // object that is derived from cppad_mixed
   bool quasi_fixed  = true;
   bool bool_sparsity = true;
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, data
   );
   mixed_object.initialize(fixed_vec, random_vec);

   // Evaluate Hessian of random likelihood w.r.t random effects
   d_sparse_rcv hes_random_obj_rcv =
      mixed_object.hes_random_obj(fixed_vec, random_vec);

   // check the Hessian
   ok &= hes_random_obj_rcv.nnz() == n_data;

   for(size_t k = 0; k < n_data; ++k)
   {  size_t r     = hes_random_obj_rcv.row()[k];
      size_t c     = hes_random_obj_rcv.col()[k];
      double v     = hes_random_obj_rcv.val()[k];
      ok          &= r == c;
      double sigma = fixed_vec[r] + double(r) / double(n_data);
      double check = 1.0 / (sigma * sigma);
      ok          &= fabs( check / v - 1.0 ) < eps;
   }
   return ok;
}
// END C++
