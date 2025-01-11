// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin sample_random}
{xrst_spell
  msg
  rng
  uu
}

Simulation the Posterior Distribution for Random Effects
########################################################

Syntax
******

| *error_msg* = *mixed_object* . ``sample_random`` (
| |tab| *sample* ,
| |tab| *fixed_vec* ,
| |tab| *random_ipopt_options* ,
| |tab| *random_lower* ,
| |tab| *random_upper* ,
| |tab| *random_in*
| )

See Also
********
:ref:`sample_fixed-name`

Prototype
*********
{xrst_literal
   // BEGIN PROTOTYPE
   // END PROTOTYPE
}

Purpose
*******
This routine draws samples from
the asymptotic posterior distribution for the
random effects given the model, the data, and the fixed effects; see
:ref:`theory@Sparse Observed Information` .

manage_gsl_rng
**************
It is assumed that
:ref:`manage_gsl_rng@get_gsl_rng` will return
a pointer to a GSL random number generator.

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

sample
******
This argument has prototype

   ``CppAD::vector<double>&`` *sample*

and its size is a multiple of
:ref:`derived_ctor@n_random` .
The input value of its elements does not matter.
We define

   *n_sample* = *sample_size* / *n_random*

If *error_msg* is empty, upon return
for ``i`` = 0 , ..., ``n_sample`` *-1* ,
``j`` = 0 , ..., ``n_random`` *-1* ,

   *sample* [ *i* * *n_random* + *j*  ]

is the *j*-th component of the *i*-th sample of the
optimal random effects.
The statistics of these samples is specified under
:ref:`sample_random@Covariance` below.

random_ipopt_options
********************
This argument has prototype

   ``const std::string&`` *random_ipopt_options*

and is the :ref:`ipopt_options-name` for optimizing the random effects.

fixed_vec
*********
This argument specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta`.

random_lower
************
This argument must have size equal to
:ref:`derived_ctor@n_random` and
specifies the lower limits for the optimization of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
The value minus infinity can be used to specify no lower limit.

random_upper
************
This argument must have size equal to
:ref:`derived_ctor@n_random` and
specifies the upper limits for the optimization of the random effect.
The value plus infinity can be used to specify no lower limit.

random_in
*********
This argument must have size equal to
:ref:`derived_ctor@n_random` and
specifies the initial value used for the optimization of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
It must hold that

   *random_lower* [ *i* ] <= *random_in* [ *i* ] <= *random_upper* [ *i* ]

for each valid index *i* .

Covariance
**********
Each sample of the random effects is an independent normal.
The mean for this distribution is the
:ref:`optimal random effects<theory@Optimal Random Effects, u^(theta)>`
:math:`\hat{u} ( \theta )`.
The variance of this distribution
is the inverse of the observed information
matrix; i.e.

.. math::

   f_{uu} [ \theta , \hat{u} ( \theta ) ] ^{-1}

This normal distribution is censored to be within the limits
*random_lower* , *random_upper* .

error_msg
*********
If *error_msg* is empty (non-empty),
:ref:`sample_random@sample`
values have been calculated (have not been calculated).
If *error_msg* is non-empty,
it is a message describing the problem.
{xrst_toc_hidden
   example/user/sample_random.cpp
}
Example
*******
The file :ref:`sample_random.cpp-name` is an example
and test of ``sample_random`` .

{xrst_end sample_random}
-----------------------------------------------------------------------------
*/

# include <Eigen/Core>
# include <Eigen/Cholesky>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <gsl/gsl_randist.h>
# include <cppad/mixed/exception.hpp>

std::string cppad_mixed::try_sample_random(
   d_vector&          sample               ,
   const std::string& random_ipopt_options ,
   const d_vector&    fixed_vec            ,
   const d_vector&    random_lower         ,
   const d_vector&    random_upper         ,
   const d_vector&    random_in            )
{  // case where there is nothing to do
   if( n_random_ == 0 )
      return "";
   //
   assert( sample.size() % n_random_ == 0   );
   assert( fixed_vec.size()    == n_fixed_  );
   assert( random_lower.size() == n_random_ );
   assert( random_upper.size() == n_random_ );
   assert( random_in.size()    == n_random_ );
   //
   // number of samples
   size_t n_sample = sample.size() / n_random_;
   //
   // optimal random effects
   d_vector random_opt;
   random_opt = optimize_random(
      random_ipopt_options, fixed_vec, random_lower, random_upper, random_in
   );
   // update the Cholesky factor corresponding to f_uu (theta, u)
   update_factor(fixed_vec, random_opt);
   //
   for(size_t i_sample = 0; i_sample < n_sample; i_sample++)
   {  // simulate a normal with mean zero and variance one
      d_vector w(n_random_);
      for(size_t j = 0; j < n_random_; j++)
         w[j] = gsl_ran_gaussian(CppAD::mixed::get_gsl_rng(), 1.0);
      //
      // set v to cholesky factor of f_uu(theta, u)^{-1} times w
      d_vector v(n_random_);
      bool ok = ldlt_ran_hes_.sim_cov(w, v);
      if( ! ok )
      {  std::string msg = "sample_random: Hessian w.r.t random effects"
            " is not positive definite";
         return msg;
      }
      //
      // add random_opt an truncate to random limits
      for(size_t j = 0; j < n_random_; j++)
      {  double samp = random_opt[j] + v[j];
         samp = std::min(samp, random_upper[j]);
         samp = std::max(samp, random_lower[j]);
         sample[i_sample * n_random_ + j] = samp;
      }
   }
   return "";
}
// --------------------------------------------------------------------------
// BEGIN PROTOTYPE
std::string cppad_mixed::sample_random(
   d_vector&          sample               ,
   const std::string& random_ipopt_options ,
   const d_vector&    fixed_vec            ,
   const d_vector&    random_lower         ,
   const d_vector&    random_upper         ,
   const d_vector&    random_in            )
// END PROTOTYPE
{  std::string error_msg = "";
   try
   {  error_msg = try_sample_random(
         sample                ,
         random_ipopt_options  ,
         fixed_vec             ,
         random_lower          ,
         random_upper          ,
         random_in
      );
   }
   catch(const std::exception& e)
   {  error_msg = "sample_random: std::exception: ";
      error_msg += e.what();
   }
   catch(const CppAD::mixed::exception& e)
   {  error_msg = e.message("sample_random");
   }
   return error_msg;
}
