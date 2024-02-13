// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin capture_xam.cpp}
{xrst_spell
   alloc
   andrew
   boolean
   logit
   ndebug
   poisson
   preprocessor
   rng
   royle
}

A Capture Example and Speed Test
################################

Syntax
******

| ``build/speed/capture_xam`` \\
| |tab| *random_seed*  \\
| |tab| *number_random*   \\
| |tab| *quasi_fixed*  \\
| |tab| *trace_optimize_fixed*  \\
| |tab| *ipopt_solve*  \\
| |tab| *bool_sparsity*  \\
| |tab| *hold_memory*  \\
| |tab| *derivative_test*  \\
| |tab| *start_near_solution*  \\
| |tab| *number_fixed_samples*  \\
| |tab| *number_locations*  \\
| |tab| *max_population*  \\
| |tab| *mean_population*  \\
| |tab| *mean_logit_probability*  \\
| |tab| *std_logit_probability*  \\
| |tab| *random_constraint*

Reference
*********
J. Andrew Royle,
Biometrics 60, 108-115 March 2004,
*N-Mixture Models for Estimating Population Size*
*from Spatially Replicated Counts.*

Command Arguments
*****************

random_seed
===========
This is a non-negative integer equal to the
seed for the random number generator,
to be specific,
:ref:`manage_gsl_rng@new_gsl_rng@s_in` used during the call to
``new_gsl_rng`` .

number_random
=============
This is a positive integer equal to the number of random effects.
This is also equal to the
number of times at which the measurements are made and is denoted by
:math:`T` below.

quasi_fixed
===========
This is either ``yes`` or ``no`` and is the value of
:ref:`derived_ctor@quasi_fixed` in the
``cppad_mixed`` derived class constructor.
The amount of memory used by the
:ref:`derived_ctor@mixed_derived` object,
after the information matrix is computed,
will be similar to after the initialization when *quasi_fixed* is no.

trace_optimize_fixed
====================
This is either ``yes`` or ``no`` .
If it is yes, a *print_level*  = 5
:ref:`trace<ipopt_trace-name>` of the fixed effects optimization
is included in the program output.
Otherwise the ipopt *print_level* is zero and
no such trace is printed.

ipopt_solve
===========
This is either ``yes`` or ``no`` .
If it is yes, the ``CppAD::ipopt::solve``
routine is used for optimizing the random effects,
otherwise ``CppAD::mixed::ipopt_random`` is used; see
:ref:`optimize_random@options@evaluation_method` .

bool_sparsity
=============
This is either ``yes`` or ``no`` .
If it is yes, boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.

hold_memory
===========
The CppAD memory allocator has a hold memory option will be set by

   ``CppAD::thread_alloc::hold_memory`` ( *hold_memory* );

where *hold_memory* is either ``yes`` or ``no`` .

derivative_test
===============
This is either ``yes`` or ``no`` .
If it is yes, the derivatives of functions used in the optimization
of the fixed effects are checked for correctness.
(This requires extra time).

start_near_solution
===================
This is either ``yes`` or ``no`` .
If it is yes, the initial point for the optimization
is the value of the fixed effects used to simulate the data.
Otherwise, the initial point is significantly different from this value.

number_fixed_samples
====================
This is a positive integer equal to the number of samples simulated
from the posterior distribution for the fixed effects using
:ref:`sample_fixed-name` .
The samples are used to approximation the
standard deviation for the optimal fixed effects.
One should use a large number of samples (at least 100)
for this purpose.

number_locations
================
This is a positive integer equal to the
number of locations at which the measurements are made; i.e.
:math:`R` in the reference.
Increasing this value increases the amount for each function evaluation,
but does not change the number of fixed or random effects.

max_population
==============
This is a positive integer equal to the
maximum value in the finite summation with respect to
population size; i.e.,
:math:`K` in the reference.
This must be greater than any of the simulated number of captures
at any location and time; i.e., and :math:`y_{i,t}`.
Also note that *max_population* does not affect the simulated
data :math:`y_{i,t}`.
A suggested value is five times *mean_population* .
The value of *max_population* is large enough if increasing it
takes more time but does not make a difference in the optimal fixed effects
(use the same value for the other arguments and same actual seed).

mean_population
===============
This is a positive floating point value equal to the
mean of the Poisson distribution for the population
used to simulate data values;
:math:`\lambda` in reference.

mean_logit_probability
======================
This is a positive floating point value equal to the
mean of the logit of the capture probability
(independent of the random effects)
used to simulate data values.
The is also equal to the mean of the random effects.

std_logit_probability
=====================
This is a positive floating point value equal to the
standard deviation of the logit of the capture probability
(independent of the random effects)
used to simulate data values.
The is also equal to the standard deviation of the random effects.

random_constraint
=================
This is either ``yes`` or ``no`` .
If it is ``no`` , there is no
:ref:`random constraint<problem@Maximum Likelihood@Random Constraints>`
for this example.
If it is ``yes`` ,
the random constraint is

.. math::

   0 = \hat{u}_0 ( \theta ) + \cdots + \hat{u}_{T-1} ( \theta )

where :math:`\hat{u} ( \theta )` is the
:ref:`optimal random effects<problem@Notation@Optimal Random Effects, u^(theta)>` .
The corresponding
:ref:`random constraint matrix<problem@Notation@Random Constraint Matrix, A>`
:math:`A` is the row vector of size :math:`T` with all ones; i.e.,

.. math::

   A = [ 1 , \cdots , 1 ] \in \B{R}^{1 \times T}

Output
******
Each output name, value pair is written in as *name* = *value*
where the amount of spaces surrounding the equal sign is not specified.
All of the pairs listed above are output.
In addition, the following name value pairs are also output.

cppad_mixed_version
===================
The ``cppad_mixed`` version number.

ldlt_cholmod
============
is the ``bin/run_cmake.sh`` configuration option
:ref:`run_cmake.sh@ldlt_cholmod` .

optimize_cppad_function
=======================
is the ``bin/run_cmake.sh`` configuration option
:ref:`run_cmake.sh@optimize_cppad_function` .

ndebug_defined
==============
is the ``NDEBUG`` preprocessor symbol defined.
This should be yes (no) if the ``bin/run_cmake.sh`` configuration option
:ref:`run_cmake.sh@build_type` is ``release``
(``debug`` ).

actual_seed
===========
If *random_seed* is zero,
the system clock, instead of *random_seed* ,
is used to seed the random number generator.
The actual random seed *actual_seed* is printed
so that you can reproduce results when *random_seed* is zero.

initialize_bytes
================
Is the amount of heap memory, in bytes,
added to the program during its :ref:`initialize-name` call.
Note that more temporary memory may have been used during this call.
In addition, only memory allocated using ``CppAD::thread_alloc`` is
included.

initialize_seconds
==================
Is the number of seconds used by the derived class :ref:`initialize-name` call.

optimize_fixed_seconds
======================
Is the number of seconds used by the call to
:ref:`optimize_fixed-name` that is used to compute the
optimal fixed effects.

optimize_random_seconds
=======================
Is the number of seconds used by a single call to
:ref:`optimize_random-name` that is used to compute the
optimal random effects.

information_mat_seconds
=======================
Is the number of seconds used by the call to
:ref:`information_mat-name` that computes the observed information matrix.

sample_fixed_seconds
====================
Is the number of seconds used by the call to
:ref:`sample_fixed-name` that computes the
:ref:`number_sample_fixed<capture_xam.cpp@Command Arguments@number_fixed_samples>`
samples for the fixed effects.

final_bytes
===========
Is final amount of heap memory, in bytes, added and retained by the program.
Only memory allocated using ``CppAD::thread_alloc`` is included.

sum_random_effects
==================
Is the sum of the optimal random effects.

mean_population_estimate
========================
Is the estimate for the
:ref:`capture_xam.cpp@Command Arguments@mean_population`
computed by :ref:`optimize_fixed-name` .

mean_logit_probability_estimate
===============================
Is the optimal estimate for the
:ref:`capture_xam.cpp@Command Arguments@mean_logit_probability`
(computed by :ref:`optimize_fixed-name` ).

std_logit_probability_estimate
==============================
Is the optimal estimate for the
:ref:`capture_xam.cpp@Command Arguments@std_logit_probability` .

mean_population_std
===================
Is the sample standard deviation of *mean_population_estimate*
(corresponding to the sample computed by :ref:`sample_fixed-name` ).

mean_logit_probability_std
==========================
Is the sample standard deviation of *mean_logit_probability_estimate* .

std_logit_probability_std
=========================
Is the sample standard deviation of *std_logit_probability_estimate* .

mean_population_ratio
=====================

| ( *mean_population_estimate* ``-`` *mean_population*  ) /
| |tab| *mean_population_std*

mean_logit_probability_ratio
============================

| ( *mean_logit_probability_estimate* ``-`` *mean_logit_probability*  ) /
| |tab| *mean_logit_probability_std*

std_probability_ratio
=====================

| ( *std_probability_estimate* ``-`` *std_probability*  ) /
| |tab| *std_probability_std*

capture_xam_ok
==============
The following conditions are checked. If they are all true,
*capture_xam_ok* is yes. Otherwise it is no.

#. *sum_random_effects* < 1 ``e-8`` || ( *random_constraint* == ``no`` )
#. *mean_population_ratio*  < 5.0
#. *mean_logit_probability_ratio*  < 5.0
#. *std_logit_probability_ratio*  < 5.0

If *capture_xam_ok* is yes, the program return value is
``0`` (no error condition).
Otherwise it is ``1`` (error condition).
{xrst_toc_hidden
   bin/capture_xam.sh
}
Example
*******
The file :ref:`capture_xam.sh-name` is an example using this program.

Notation
********

.. list-table::
   :widths: auto

   * - :math:`R`
     - number of sampling locations; i.e., *number_locations* .
   * - :math:`T`
     - number of sampling times; i.e., *number_random* .
   * - :math:`K`
     - maximum population in truncation of infinite summation; i.e.,
       *max_population*
   * - :math:`\theta_0`
     - mean for :math:`N_i` given :math:`\theta`
       (:math:`\lambda` in reference); i.e.,
       *mean_population* .
   * - :math:`\theta_1`
     - mean of logit of capture probability; i.e.,
       *mean_logit_probability* .
   * - :math:`\theta_2`
     - standard deviation of logit of capture probability; i.e.,
       *std_logit_probability* .
   * - :math:`N_i`
     - size of the population at *i*-th location
   * - :math:`y_{i,t}`
     - number of captures at location :math:`i` and time :math:`t`
       (:math:`n_{i,t}` in reference)
   * - :math:`y_i`
     - is the vector of captures at location :math:`i`
       :math:`( y_{i,0} , \ldots , y_{i, T-1} )`.
   * - :math:`M_i`
     - maximum of captures at *i*-th location
       (:math:`\max_t n_{i,t}` in reference)
   * - :math:`q_t`
     - capture probability at :math:`t` same for all locations
       (:math:`p_{i,t}` in reference)
   * - :math:`u_t`
     - random effect for each sampling time

p(y_it|N_i,q_t)
***************
We use a binomial distribution to model the
probability of :math:`y_{i,t}` given :math:`N_i` and :math:`q_t`; i.e,

.. math::

   \B{p} ( y_{i,t} | N_i , q_t )
   =
   \left( \begin{array}{c} N_i \\ y_{i,t} \end{array} \right)
   q_t^{y(i,t)} \left( 1 - q_t \right)^{y(i,t)}

Furthermore, we assume that this probability
is independent for each :math:`(i, t)`.

p(N_i|theta)
************
We use a Poisson distribution to model the
probability of :math:`N_i` given :math:`\theta_0`; i.e.,

.. math::

   \B{p} ( N_i | \theta  )
   =
   \theta_0^{N(i)} \frac{ \exp[ - \theta_0 ] }{ N_i ! }

We assume these this probability
is independent for each :math:`i`.

q_t(theta,u)
************
Section 2.4 of the
:ref:`capture_xam.cpp@Reference` suggests a
covariate model for the probability of capture.
We use a similar model defined by

.. math::

   \R{logit} ( q_t ) = u_t + \theta_1

It follows that

.. math::

   q_t( \theta , u) = [ 1 + \exp(- u_t - \theta_1 ) ]^{-1}

M_i
***
We define the vector of maximum measurement for each location by

.. math::

   M_i = \max \left\{ y_{i,0} , \cdots , y_{i, T-1} \right\}

p(y_i|theta,u)
**************
The probability for :math:`y_i` (the captures at location :math:`i`)
given :math:`N_i`, :math:`\theta`, and :math:`u` is

.. math::

   \B{p}( y_i | N_i, \theta , u )
   =
   \prod_{t=0}^{T-1}
   \left( \begin{array}{c} {N(i)} \\ y_{i,t} \end{array} \right)
      q_t ( \theta , u)^{y(i,t)}
      \left( 1 - q_t( \theta , u) \right)^{y(i,t)}

We do not know the population at each location :math:`N_i`,
but instead have a Poisson prior for :math:`N_i`.
We sum with respect to the possible values for :math:`N_i`
to get the probability of :math:`y_i` given
:math:`\theta` and :math:`u`.

.. math::

   \B{p}( y_i | \theta , u )
   =
   \sum_{k=0}^K \B{p}( y_i | N_i=k, \theta , u ) \B{p}( N_i=k | \theta )

where :math:`k` is the possible values for :math:`N_i`.
Note that :math:`K` should be plus infinity, but we use a fixed
finite value for :math:`K` as an approximation for the infinite sum.

p(y|theta,u)
************
Our model for the likelihood of the data at all the locations,
given the fixed and random effects, is

.. math::

   \B{p}( y | \theta , u )
   =
   \prod_{i=0}^{R-1} \B{p}( y_i | \theta , u )

Expressed in terms of fixed effects :math:`\theta`,
the random effects :math:`u`,
and the data :math:`y`, this is

.. math::

   \B{p}( y | \theta , u )
   =
   \prod_{i=0}^{R-1}
   \left[
   \sum_{k=0}^{K-1}
   \theta_0^k \frac{ \exp[ - \theta_0 ] }{ k ! }
   \prod_{t=0}^{T-1}
   \left( \begin{array}{c} {k} \\ y_{i,t} \end{array} \right)
      q_t ( \theta , u)^{y(i,t)}
      \left( 1 - q_t( \theta , u) \right)^{y(i,t)}
   \right]

In ``cppad_mixed`` notation, this specifies the
:ref:`random data density<problem@Notation@Random Data Density, p(y|theta,u)>` .

p(u|theta)
**********
We use a normal distribution, with mean zero and standard deviation
:math:`\theta_2`,
for the distribution of the random effects :math:`u`
given the fixed effects :math:`\theta`; i.e.,

.. math::

   \B{p} ( u | \theta )
   =
   \prod_{t=0}^{T-1}
      \frac{1}{ \theta_2 \sqrt{ 2 \pi } }
         \exp \left[ - \frac{1}{2} \frac{ u_t^2 }{ \theta_2^2 } \right]

In ``cppad_mixed`` notation, this specifies the
:ref:`random prior density<problem@Notation@Random Prior Density, p(u|theta)>` .

p(theta)
********
For this example there is no
:ref:`fixed prior density<problem@Notation@Fixed Prior Density, p(theta)>`
:math:`\B{p}(\theta)`.

p(z|theta)
**********
For this example there is no
:ref:`fixed data density<problem@Notation@Fixed Data Density, p(z|theta)>`
:math:`\B{p}(z | \theta)`.

c(theta)
********
For this example there is no
:ref:`fixed constraint function<problem@Notation@Fixed Constraint Function, c(theta)>`
:math:`c( \theta )`.

Source Code
***********
{xrst_literal
   // BEGIN C++
   // END C++
}

{xrst_end capture_xam.cpp}
-----------------------------------------------------------------------------
*/
// BEGIN C++
# include <ctime>
# include <gsl/gsl_randist.h>
# include <cppad/utility/vector.hpp>
# include <cppad/utility/elapsed_seconds.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <cppad/mixed/configure.hpp>

namespace { // BEGIN_EMPTY_NAMESPACE
//
using std::exp;
using std::log;
//
using CppAD::mixed::s_vector;
using CppAD::mixed::d_vector;
using CppAD::mixed::d_sparse_rcv;
//
// Convert size_t to string adding commas every three digits
std::string size_t2string(size_t value )
{  std::string raw_string = CppAD::to_string(value);
   std::string result = "";
   size_t n_raw = raw_string.size();
   for(size_t i = 0; i < n_raw; i++)
   {  if( i != 0 && (n_raw - i) % 3 == 0 )
         result += ",";
      result += raw_string[i];
   }
   return result;
}

// simulate data, y
void simulate(
   bool                   random_constraint ,
   size_t                 R                 ,
   size_t                 T                 ,
   const d_vector&        theta             ,
   s_vector&              y                 )
{  assert( theta.size() == 3 );
   assert( y.size() == R * T );
   //
   // random number generator
   gsl_rng* rng = CppAD::mixed::get_gsl_rng();
   //
   // simulate population sizes
   d_vector N(R);
   double mu =  theta[0];
   for(size_t i = 0; i < R; i++)
      N[i] = gsl_ran_poisson(rng, mu );
   //
   // simulate random effects
   d_vector u(T);
   double sigma = theta[2];
   double sum   = 0.0;
   for(size_t t = 0; t < T; t++)
   {  u[t] = gsl_ran_gaussian(rng, sigma);
      sum += u[t];
   }
   // adjust the random effects when using random constraint
   if( random_constraint )
   {  for(size_t t = 0; t < T; t++)
      {  // simulate mean zero values
         u[t] = u[t] - sum / double(T);
         //
         // correct for loss of one degree of freedom
         u[t] = u[t] * sqrt( double(T) / double(T-1) );
      }
   }
   // simulate data
   for(size_t i = 0; i < R; i++)
   {  for(size_t t = 0; t < T; t++)
      {  // probability of capture
         double ex = exp( - u[t] - theta[1] );
         double q = 1.0 /( 1.0  + ex );
         y[ i * T + t ] = gsl_ran_binomial(rng, q, (unsigned int)N[i] );
      }
   }
   //
   return;
}

// cppad_mixed derived class
class mixed_derived : public cppad_mixed {
private:
   const size_t          R_; // number of locations
   const size_t          T_; // number of times
   size_t                K_; // max used when summing over population size
   const s_vector&       y_; // reference to data values
   // -----------------------------------------------------------------
   // set by constructor and then effectively const
   s_vector        M_;       // max number of captures at each location
   d_vector        logfac_;  // logfac_[k] = log( k! )
   d_vector        log_pik_; // used by implement_ran_likelihood
// ------------------------------------------------------------------------
public:
   // constructor
   mixed_derived(
      size_t                  R             ,
      size_t                  T             ,
      size_t                  K             ,
      bool                    quasi_fixed   ,
      bool                    bool_sparsity ,
      const d_sparse_rcv&     A_rcv         ,
      s_vector&               y             ,
      d_vector&               fixed_in      ,
      d_vector&               random_in     ) :
      // n_fixed = 3, n_random = T
      cppad_mixed(
         3, T, quasi_fixed, bool_sparsity, A_rcv
      )                ,
      R_(R)            ,
      T_(T)            ,
      K_(K)            ,
      y_(y)
   {  assert( fixed_in.size() == 3 );
      assert( random_in.size() == T );
      //
      // set M_
      M_.resize(R);
      for(size_t i = 0; i < R; i++)
      {  M_[i] = 0;
         for(size_t t = 0; t < T; t++)
            M_[i] = std::max( M_[i], y[ i * T + t] );
         assert( M_[i] < K_ );
      }
      //
      // set logfac_
      logfac_.resize(K_);
      logfac_[0]  = 0.0;
      logfac_[1]  = 0.0;
      for(size_t k = 2; k < K_; k++)
         logfac_[k] = log( double(k) ) + logfac_[k-1];
      //
      // set log_pik_
      log_pik_.resize(R_ * K_);
      for(size_t ell = 0; ell < R_ * K_; ell++)
         log_pik_[ell] = 0.0;
      template_ran_likelihood(fixed_in, random_in, log_pik_);
   }
   // implementaion of ran_likelihood, used with Float = double and a1_double
   template <class Float>
   CppAD::vector<Float> template_ran_likelihood(
      const CppAD::vector<Float>&  theta   ,
      const CppAD::vector<Float>&  u       ,
      CppAD::vector<Float>&        log_pik )
   {  using CppAD::vector;
      //
      assert( log_pik.size() == R_* K_ );
      vector<Float> vec(1);
      //
      Float one( 1.0 );
      Float two( 2.0 );
      Float pi2( 8.0 * std::atan(1.0) );
      Float eps = Float( 10.0 * std::numeric_limits<double>::epsilon() );
      //  ------------------------------------------------------------
      // log [ p(u | theta) ]
      //  ------------------------------------------------------------
      Float sig = theta[2];
      vec[0] -= Float(T_) * ( two * log(sig) + log( pi2 ) );
      for(size_t t = 0; t < T_; t++)
      {  Float w = u[t] / ( sig + eps );
         vec[0] -= w * w;
      }
      vec[0] /= two;
      //  ------------------------------------------------------------
      // log [ p(y | theta, u) ]
      //  ------------------------------------------------------------
      //
      // y_{i,t} * log( qt )
      vector<Float> yit_log_q(R_ * T_);
      // log( 1.0 - qt )
      vector<Float> log_1q(T_);
      for(size_t t = 0; t < T_; t++)
      {  Float ex  = exp( - u[t] - theta[1] );
         Float q   = one / (one + ex );
         log_1q[t] = log(one - q + eps);
         for(size_t i = 0; i < R_; i++)
            yit_log_q[ i * T_ + t ]  = Float(y_[i*T_+t]) * log(q + eps);
      }
      //
      // log( theta[0]^k * exp( - theta[0] ) / k! )
      vector<Float> log_poisson(K_);
      for(size_t k = 0; k < K_; k++)
      {  log_poisson[k] = log( theta[0] + eps ) * Float(k);
         log_poisson[k] -= theta[0];
         log_poisson[k] -= logfac_[k];
      }
      //
      // log[ p(y|theta, u) ] to vec[0]
      for(size_t i = 0; i < R_; i++)
      {  // compute maximum of log_pik for this i
         Float max_log_pik = log_pik[ i * K_ + M_[i] ];
         for(size_t k = M_[i] + 1; k < K_; k++)
         {  if( log_pik[ i * K_ + k ] > max_log_pik )
               max_log_pik = log_pik[ i * K_ + k ];
         }
         //
         // initialize p(y_i|theta,u)
         Float p_i = Float(0.0);
         for(size_t k = M_[i]; k < K_; k++)
         {  // compute p(y_i|N_i=k,theta,u) p(N_i=k|theta)
            //
            // initialize sum that needs to be calculated with Float
            // p(N_i=k|theta)
            Float float_sum = log_poisson[k];
            //
            // initialize sum that does not need to use Float
            double double_sum = 0.0;
            //
            // now compute terms that depend on t
            for(size_t t = 0; t < T_; t++)
            {  size_t yit = y_[ i * T_ + t ];
               //
               // k choose yit = k! / [ yit! * (k - yit)! ]
               // log [ (k choose yit) ]
               double_sum += logfac_[k] - logfac_[yit] - logfac_[k - yit];
               //
               // log ( qt^yit )
               float_sum += yit_log_q[i * T_ + t];
               //
               // log [ (1 - qt)^(k - yit) ]
               float_sum += Float(k - yit) * log_1q[t];
            }
            // log[ p(y_i|N_i=k,theta,u) p(N_i=k|theta)
            // (output elements of log_pik are not affected by its input)
            log_pik[ i * K_ + k ] = float_sum + double_sum;
            //
            // normalization avoids exponential overflow
            // p_i += p(y_i|N_i=k,theta,u) p(N_i=k|theta) / exp(max_log_pik)
            p_i += exp( log_pik[ i * K_ + k ] - max_log_pik );
         }
         // vec[0] += log[ p(y_i|theta,u) ]
         vec[0] += log( p_i );
         // following term does not depend on the fixed or random effects
         // vec[0] += Float(K_ - M_[i]) * max_log_pik;
      }
      // - log [ p(y|theta,u) p(u|theta) ]
      vec[0] = - vec[0];
      //
      // result may be inifite or nan when only computing log_pik
      // (initial log_pik is used to compute normalization factors)
      return vec;
   }
// ------------------------------------------------------------------------
public:
   // a1_vector ran_likelihood
   virtual a1_vector ran_likelihood(
      const a1_vector& fixed_vec  ,
      const a1_vector& random_vec )
   {  a1_vector log_pik(R_ * K_), vec(1);
      for(size_t ell = 0; ell < R_ * K_; ell++)
         log_pik[ell] = a1_double(log_pik_[ell]);
      vec = template_ran_likelihood(fixed_vec, random_vec, log_pik);
      // make sure result is finite
# ifndef NDEBUG
      double inf = std::numeric_limits<double>::infinity();
# endif
      assert( vec[0] < + a1_double(inf) );
      assert( vec[0] > - a1_double(inf) );
      //
      return vec;
   }
};
template <class Value>
void label_print(const char* label, const Value& value)
{  std::cout << std::setw(35) << std::left << label;
   std::cout << " = " << value << std::endl;
}
void label_print(const char* label, const double& value)
{  std::cout << std::setw(35) << std::left << label;
   int n_digits = 3 + int( std::log10(value) + 1e-9 );
   std::cout << " = " << std::setprecision(n_digits) << value << std::endl;
}
void bool_print(const char* label, bool value)
{  std::cout << std::setw(35) << std::left << label;
   if( value )
      std::cout << " = yes";
   else
      std::cout << " = no";
   std::cout << std::endl;
}

} // END_EMPTY_NAMESPACE

int main(int argc, const char *argv[])
{  bool ok = true;
   using std::endl;
   using std::string;
   //
   const char* arg_name[] = {
      "random_seed",
      "number_random",
      "quasi_fixed",
      "trace_optimize_fixed",
      "ipopt_solve",
      "bool_sparsity",
      "hold_memory",
      "derivative_test",
      "start_near_solution",
      "number_fixed_samples",
      "number_locations",
      "max_population",
      "mean_population",
      "mean_logit_probability",
      "std_logit_probability",
      "random_constraint"
   };
   size_t n_arg = sizeof(arg_name)/sizeof(arg_name[0]);
   //
   if( size_t(argc) != 1 + n_arg )
   {  // print usage error message
      std::cerr << "expected " << n_arg << " arguments and found "
      << argc - 1 << "\nusage: " << argv[0];
      for(size_t i = 0; i < n_arg; i++)
         std::cerr << " \\ \n\t" << arg_name[i];
      std::cerr << "\n";
      std::exit(1);
   }
   //
   // get command line arguments
   assert( n_arg == 16 );
   size_t iarg = 1;
   size_t random_seed              = std::atoi( argv[iarg++] );
   size_t number_random            = std::atoi( argv[iarg++] );
   string quasi_fixed_str          = argv[iarg++];
   string trace_optimize_fixed_str = argv[iarg++];
   string ipopt_solve_str          = argv[iarg++];
   string bool_sparsity_str        = argv[iarg++];
   string hold_memory_str          = argv[iarg++];
   string derivative_test_str      = argv[iarg++];
   string start_near_solution_str  = argv[iarg++];
   //
   size_t number_fixed_samples   = std::atoi( argv[iarg++] );
   size_t number_locations       = std::atoi( argv[iarg++] );
   size_t max_population         = std::atoi( argv[iarg++] );
   double mean_population        = std::atof( argv[iarg++] );
   double mean_logit_probability = std::atof( argv[iarg++] );
   double std_logit_probability  = std::atof( argv[iarg++] );
   string random_constraint_str  = argv[iarg++];
   //
   bool random_constraint    =  random_constraint_str == "yes";
   bool quasi_fixed          =  quasi_fixed_str == "yes";
   bool trace_optimize_fixed =  trace_optimize_fixed_str == "yes";
   bool ipopt_solve          =  ipopt_solve_str == "yes";
   bool bool_sparsity        =  bool_sparsity_str == "yes";
   bool hold_memory          =  hold_memory_str == "yes";
   bool derivative_test      =  derivative_test_str == "yes";
   bool start_near_solution  =  start_near_solution_str == "yes";
   //
   // hold memory setting
   CppAD::thread_alloc::hold_memory(hold_memory);
   //
   // print the command line arugments with labels for each value
   for(size_t i = 0; i < n_arg; i++)
      label_print(arg_name[i], argv[1+i]);
   //
   assert( number_fixed_samples > 0 );
   assert( number_locations > 0 );
   assert( number_random > 0 );
   assert( max_population > 0 );
   assert( mean_population > 0.0 );
   assert( std_logit_probability > 0.0 );
   assert( random_constraint || random_constraint_str=="no" );
   assert( quasi_fixed || quasi_fixed_str=="no" );
   assert( trace_optimize_fixed || trace_optimize_fixed_str=="no" );
   assert( ipopt_solve || ipopt_solve_str=="no" );
   assert( bool_sparsity || bool_sparsity_str =="no" );
   assert( hold_memory || hold_memory_str =="no" );
   assert( derivative_test || derivative_test_str =="no" );
   assert( start_near_solution || start_near_solution_str =="no" );
   //
   // configuration options
   label_print("cppad_mixed_version",    CPPAD_MIXED_VERSION);
   bool_print("ldlt_cholmod",            CPPAD_MIXED_LDLT_CHOLMOD);
   bool_print("optimize_cppad_function", CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION);
# ifdef NDEBUG
   bool_print("ndebug_defined",          true);
# else
   bool_print("ndebug_defined",          false);
# endif
   //
   // Get actual seed (different when random_seed is zero).
   size_t actual_seed     = CppAD::mixed::new_gsl_rng( random_seed );
   label_print("actual_seed", actual_seed);
   //
   // number of fixed and random effects
   size_t n_fixed  = 3;
   size_t n_random = number_random;
   //
   // number of random effects
   size_t T = number_random;
   size_t R = number_locations;
   size_t K = max_population;
   //
   d_vector theta_sim(n_fixed);
   theta_sim[0] = mean_population;
   theta_sim[1] = mean_logit_probability;
   theta_sim[2] = std_logit_probability;
   //
   // simulate y
   s_vector y(R * T);
   simulate(random_constraint, R, T, theta_sim, y);

   // lower and upper limits
   d_vector fix_constraint_lower, fix_constraint_upper;
   d_vector theta_lower(n_fixed), theta_in(n_fixed), theta_upper(n_fixed);
   //
   // theta[0] is estimate of mean_population
   theta_lower[0] = 0.0;
   theta_in[0]    = mean_population / 2.0;
   theta_upper[0] = 2.0 * mean_population;
   //
   // theta[1] is estimate of mean_logit_probability
   theta_lower[1] = mean_logit_probability - 2.0;
   theta_in[1]    = mean_logit_probability - 0.5;
   theta_upper[1] = mean_logit_probability + 2.0;
   //
   // theta[2] is estimate of std_logit_probability
   theta_lower[2] = std_logit_probability / 10.;
   theta_in[2]    = std_logit_probability / 2.0;
   theta_upper[2] = std_logit_probability * 10.;

   // check if we are starting at the simulation values
   if( start_near_solution )
   {  for(size_t j = 0; j < n_fixed; j++)
         theta_in[j] = theta_sim[j];
   }

   // random constraints
   CppAD::mixed::sparse_rc A_pattern;
   if( random_constraint )
   {  A_pattern.resize(1, T, T);
      for(size_t t = 0; t < T; t++)
         A_pattern.set(t, 0, t);
   }
   d_sparse_rcv A_rcv( A_pattern );
   for(size_t t = 0; t < A_rcv.nc(); t++)
      A_rcv.set(t, 1.0);

   // initialize random effects to start optimization at
   d_vector u_in(T);
   for(size_t t = 0; t < T; t++)
      u_in[t] = 0.0;
   //
   // create derived object
   mixed_derived mixed_object(
      R, T, K, quasi_fixed, bool_sparsity, A_rcv, y, theta_in, u_in
   );
   //
   // start timing of initialize
   size_t thread         = CppAD::thread_alloc::thread_num();
   size_t start_bytes    = CppAD::thread_alloc::inuse(thread);
   double start_seconds  = CppAD::elapsed_seconds();
   //
   mixed_object.initialize(theta_in, u_in);
   //
   double end_seconds = CppAD::elapsed_seconds();
   size_t end_bytes   = CppAD::thread_alloc::inuse(thread);
   //
   // print amoumt of memory added to mixed_object during initialize
   // (use commans to separate every three digits).
   string initialize_bytes = size_t2string(end_bytes - start_bytes);
   label_print("initialize_bytes", initialize_bytes);
   label_print("initialize_seconds", end_seconds - start_seconds);
   //
   // ipopt options for optimizing the random effects
   string random_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           none\n"
      "Numeric tol                       1e-9\n"
   ;
   if( ipopt_solve ) random_ipopt_options +=
      "String evaluation_method         ipopt_solve\n";
   //
   // ipopt options for optimizing the fixd effects
   string fixed_ipopt_options =
      "String  sb                        yes\n"
      "Numeric tol                       1e-9\n"
      "Integer max_iter                  40\n"
   ;
   if( trace_optimize_fixed )
      fixed_ipopt_options += "Integer print_level         5\n";
   else
      fixed_ipopt_options += "Integer print_level         0\n";
   //
   if( derivative_test )
      fixed_ipopt_options += "String derivative_test      first-order\n";
   else
      fixed_ipopt_options += "String derivative_test      none\n";
   //
   double inf = std::numeric_limits<double>::infinity();
   d_vector u_lower(n_random), u_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  u_lower[i] = -inf;
      u_upper[i] = +inf;
   }
   //
   // optimize fixed effects
   start_seconds = CppAD::elapsed_seconds();
   if( trace_optimize_fixed )
      std::cout << endl;
   d_vector theta_scale = theta_in;
   CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
      fixed_ipopt_options,
      random_ipopt_options,
      theta_lower,
      theta_upper,
      fix_constraint_lower,
      fix_constraint_upper,
      theta_scale,
      theta_in,
      u_lower,
      u_upper,
      u_in
   );
   end_seconds = CppAD::elapsed_seconds();
   if( trace_optimize_fixed )
      std::cout << endl;
   label_print("optimize_fixed_seconds", end_seconds - start_seconds);
   //
   // estimate of fixed effects
   d_vector theta_out = solution.fixed_opt;
   //
   // estimate of random effects
   start_seconds = CppAD::elapsed_seconds();
   d_vector u_out     = mixed_object.optimize_random(
      random_ipopt_options,
      theta_out,
      u_lower,
      u_upper,
      u_in
   );
   end_seconds = CppAD::elapsed_seconds();
   label_print("optimize_random_seconds", end_seconds - start_seconds);
   //
   // information matrix
   start_seconds = CppAD::elapsed_seconds();
   d_sparse_rcv hes_fixed_obj_rcv = mixed_object.hes_fixed_obj(
      solution.fixed_opt, u_out
   );
   end_seconds = CppAD::elapsed_seconds();
   label_print("information_mat_seconds", end_seconds - start_seconds);
   //
   // sample approximate posteroior for fixed effects
   start_seconds = CppAD::elapsed_seconds();
   d_vector sample( number_fixed_samples * n_fixed );
   mixed_object.sample_fixed(
      sample,
      hes_fixed_obj_rcv,
      solution,
      theta_lower,
      theta_upper
   );
   end_seconds = CppAD::elapsed_seconds();
   label_print("sample_fixed_seconds", end_seconds - start_seconds);
   //
   end_bytes          = CppAD::thread_alloc::inuse(thread);
   string final_bytes = size_t2string(end_bytes - start_bytes);
   label_print("final_bytes", final_bytes);
   //
   // sum of random effects
   double sum_random_effects = 0.0;
   for(size_t j = 0; j < n_random; j++)
      sum_random_effects += u_out[j];
   if( random_constraint )
      ok &= std::fabs( sum_random_effects ) < 1e-8;
   label_print("sum_random_effects", sum_random_effects);
   //
   // sample_std
   d_vector sample_std(n_fixed);
   for(size_t j = 0; j < n_fixed; j++)
      sample_std[j] = 0.0;
   for(size_t i = 0; i < number_fixed_samples; i++)
   {  for(size_t j = 0; j < n_fixed; j++)
      {  double diff    = sample[i * n_fixed + j] - theta_out[j];
         sample_std[j] += diff * diff;
      }
   }
   // estimate_ratio
   d_vector estimate_ratio(n_fixed);
   for(size_t j = 0; j < n_fixed; j++)
   {  sample_std[j] = std::sqrt(sample_std[j]/double(number_fixed_samples));
      // check if results results are reasonable
      estimate_ratio[j] = ( theta_out[j] - theta_sim[j] ) / sample_std[j];
      ok  &= std::fabs(estimate_ratio[j]) < 10.0;
   }
   // estimates
   label_print("mean_population_estimate", theta_out[0]);
   label_print("mean_logit_probability_estimate", theta_out[1]);
   label_print("std_logit_probability_estimate", theta_out[2]);
   // simulation sample standard deviations
   label_print("mean_population_std", sample_std[0]);
   label_print("mean_logit_probability_std", sample_std[1]);
   label_print("std_logit_probability_std", sample_std[2]);
   // ratios, error conditon is ratio >= 5.0
   label_print("mean_population_ratio", estimate_ratio[0]);
   label_print("mean_logit_probability_ratio", estimate_ratio[1]);
   label_print("std_logit_probability_ratio", estimate_ratio[2]);
   //
   if( ok )
      label_print("capture_xam_ok", "yes");
   else
      label_print("capture_xam_ok", "no");
   //
   // make sure CppAD::thread_alloc is no longer holding memory
   CppAD::thread_alloc::hold_memory(false);
   // free memory allocated by new_gsl_rng
   CppAD::mixed::free_gsl_rng();
   if( ok )
      return 0;
   return 1;
}
// END C++
