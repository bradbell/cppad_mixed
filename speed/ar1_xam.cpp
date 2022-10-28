// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
$begin ar1_xam.cpp$$
$spell
   cmake
   bool
   ar1_xam
   gsl_rng
   ipopt
   cppad
   CppAD
   alloc
   cholesky
   ldlt_cholmod
   ndebug
$$

$section A First Order Auto-Regressive Example and Speed Test$$

$head Syntax$$
$codei%./ar1_xam \
   %random_seed% \
   %number_random% \
   %quasi_fixed% \
   %trace_optimize_fixed% \
   %ipopt_solve% \
   %bool_sparsity% \
   %hold_memory% \
   %derivative_test% \
   %start_near_solution%
%$$

$head Problem$$

$subhead Data$$
For $latex t = 0 , \ldots , T - 1$$,
$latex y_t = (1 + t) + e_t $$,
where $latex e_t \sim \B{N}( 0, \sigma_y^2 )$$.


$subhead p( y_t | u , theta )$$
For $latex t = 0 , \ldots , T - 1$$,
$latex y_t \sim \B{N}( u_t , \sigma_y^2 )$$.

$subhead p( u | theta )$$
For $latex t = 0$$, $latex u_t \sim \B{N}( 0 , \theta_0^2 )$$,
and for $latex t = 1 , \ldots , T - 1$$,
$latex u_t - u_{t-1} \sim \B{N}( 0 , \theta_0^2 )$$.


$head Command Arguments$$

$subhead random_seed$$
This is a non-negative integer equal to the
seed for the random number generator,
to be specific,
$cref/s_in/manage_gsl_rng/new_gsl_rng/s_in/$$ used during the call to
$code new_gsl_rng$$.

$subhead number_random$$
This is a positive integer specifying the number of random effects.
This is also the number of time points and number of data values.

$subhead quasi_fixed$$
This is either $code yes$$ or $code no$$ and is the value of
$cref/quasi_fixed/derived_ctor/quasi_fixed/$$ in the
$code cppad_mixed$$ derived class constructor.
The amount of memory used by the
$cref/mixed_derived/derived_ctor/mixed_derived/$$ object,
after the information matrix is computed,
will be similar to after the initialization when $icode quasi_fixed$$ is no.

$subhead trace_optimize_fixed$$
This is either $code yes$$ or $code no$$.
If it is yes, a $icode%print_level% = 5%$$
$cref/trace/ipopt_trace/$$ of the fixed effects optimization
is included in the program output.
Otherwise the ipopt $icode print_level$$ is zero and
no such trace is printed.

$subhead ipopt_solve$$
This is either $code yes$$ or $code no$$.
If it is yes, the $code CppAD::ipopt::solve$$
routine is used for optimizing the random effects,
otherwise $code CppAD::mixed::ipopt_random$$ is used; see
$cref/evaluation_method/optimize_random/options/evaluation_method/$$.

$subhead bool_sparsity$$
This is either $code yes$$ or $code no$$.
If it is yes, boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.

$subhead hold_memory$$
The CppAD memory allocator has a hold memory option will be set by
$codei%
   CppAD::thread_alloc::hold_memory(%hold_memory%);
%$$
where $icode hold_memory$$ is either $code yes$$ or $code no$$.

$subhead derivative_test$$
This is either $code yes$$ or $code no$$.
If it is yes, the derivatives of functions used in the optimization
of the fixed effects are checked for correctness.
(This requires extra time).

$subhead start_near_solution$$
This is either $code yes$$ or $code no$$.
If it is yes, the initial point for the optimization
is the value of the fixed effects used to simulate the data.
Otherwise, the initial point is significantly different from this value.


$head Output$$
Each output name, value pair is written in as $icode%name% = %value%$$
where the amount of spaces surrounding the equal sign is not specified.
All of the pairs listed above are output.
In addition, the following name value pairs are also output.

$subhead cppad_mixed_version$$
The $code cppad_mixed$$ version number.

$subhead ldlt_cholmod$$
is the $code bin/run_cmake.sh$$ configuration option
$cref/ldlt_cholmod/run_cmake.sh/ldlt_cholmod/$$.

$subhead optimize_cppad_function$$
is the $code bin/run_cmake.sh$$ configuration option
$cref/optimize_cppad_function/run_cmake.sh/optimize_cppad_function/$$.

$subhead ndebug_defined$$
is the $code NDEBUG$$ preprocessor symbol defined.
This should be yes (no) if the $code bin/run_cmake.sh$$ configuration option
$cref/build_type/run_cmake.sh/build_type/$$ is $code release$$
($code debug$$).

$subhead actual_seed$$
If $icode random_seed$$ is zero,
the system clock, instead of $icode random_seed$$,
is used to seed the random number generator.
The actual random seed $icode actual_seed$$ is printed
so that you can reproduce results when $icode random_seed$$ is zero.

$subhead initialize_bytes$$
Is the amount of heap memory, in bytes,
added to the program during its $cref initialize$$ call.
Note that more temporary memory may have been used during this call.
In addition, only memory allocated using $code CppAD::thread_alloc$$ is
included.

$subhead initialize_seconds$$
Is the number of seconds used by the derived class $cref initialize$$ call.

$subhead optimize_fixed_seconds$$
Is the number of seconds used by the call to
$cref optimize_fixed$$ that is used to compute the
optimal fixed effects.

$subhead optimize_random_seconds$$
Is the number of seconds used by a single call to
$cref optimize_random$$ that is used to compute the
optimal random effects.

$subhead information_mat_seconds$$
Is the number of seconds used by the call to
$cref information_mat$$ that computes the observed information matrix.

$subhead sample_fixed_seconds$$
Is the number of seconds used by the call to
$cref sample_fixed$$ that computes the
$cref/number_sample_fixed
   /capture_xam.cpp/Command Arguments/number_fixed_samples/$$
samples for the fixed effects.

$subhead final_bytes$$
Is final amount of heap memory, in bytes, added and retained by the program.
Only memory allocated using $code CppAD::thread_alloc$$ is included.

$subhead theta_0_estimate$$
Is the optimal estimate for $latex \theta_0$$; see the
$cref/problem/ar1_xam.cpp/Problem/$$ definition.

$subhead ar1_xam_ok$$
If this program passes it's correctness test,
$icode ar1_xam_ok$$ is yes and the program return code is $code 0$$.
Otherwise $icode ar1_xam_ok$$ it is no and the return code is $code 1$$.

$children%bin/ar1_xam.sh
%$$
$head Example$$
The file $cref ar1_xam.sh$$ is an example using this program.

$head Source Code$$
$code
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$


$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <gsl/gsl_randist.h>
# include <Eigen/Dense>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/manage_gsl_rng.hpp>
# include <cppad/mixed/configure.hpp>

namespace {
   using std::cout;
   using std::string;
   using CppAD::log;
   using CppAD::AD;
   //
   using CppAD::mixed::d_sparse_rcv;
   using CppAD::mixed::a1_double;
   using CppAD::mixed::d_vector;
   //
   // The constant sigma_y
   const double sigma_y = 0.1;
   //
   class mixed_derived : public cppad_mixed {
   private:
      const size_t          n_fixed_;
      const size_t          n_random_;
      const d_vector&       y_;
   // ----------------------------------------------------------------------
   public:
      // constructor
      mixed_derived(
         size_t                 n_fixed       ,
         size_t                 n_random      ,
         bool                   quasi_fixed   ,
         bool                   bool_sparsity ,
         const d_sparse_rcv&    A_rcv         ,
         const d_vector&        y             ) :
         cppad_mixed(
            n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
         )                     ,
         n_fixed_(n_fixed)     ,
         n_random_(n_random)   ,
         y_(y)
      {  assert( n_fixed == 1);
         assert( y_.size() == n_random_ );
      }
      // -------------------------------------------------------------------
      // implementation of ran_likelihood
      // Note that theta[2] is not used
      template <typename Vector>
      Vector template_ran_likelihood(
         const Vector& theta  ,
         const Vector& u      )
      {  typedef typename Vector::value_type scalar;

         assert( theta.size() == n_fixed_ );
         assert( u.size() == y_.size() );
         Vector vec(1);
         //
         // initialize summation
         vec[0] = scalar(0.0);
         //
         // sqrt_2pi = CppAD::sqrt(8.0 * CppAD::atan(1.0) );
         //
         for(size_t i = 0; i < n_random_; i++)
         {  scalar sigma  = scalar(sigma_y);
            scalar res    = (y_[i] - u[i]) / sigma;
            //
            // p(y_i | u, theta)
            vec[0] += res * res / scalar(2.0);
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);
            //
            // p(u_i | theta)
            sigma = theta[0];
            if( i == 0 )
               res = u[i] / sigma;
            else
               res = ( u[i] - u[i-1] ) / sigma;
            vec[0] += log(sigma) + res * res / scalar(2.0);
            // following term does not depend on fixed or random effects
            // vec[0] += log(sqrt_2pi);
         }
         return vec;
      }
      // a1_vector version of ran_likelihood
      virtual a1_vector ran_likelihood(
         const a1_vector& fixed_vec, const a1_vector& random_vec
      )
      {  return template_ran_likelihood( fixed_vec, random_vec ); }
   };
   template <class Value>
   void label_print(const char* label, const Value& value)
   {  cout << std::setw(35) << std::left << label;
      cout << " = " << value << std::endl;
   }
   void label_print(const char* label, const double& value)
   {  cout << std::setw(35) << std::left << label;
      int n_digits = 5 + int( std::log10(value) + 1e-9 );
      cout << " = " << std::setprecision(n_digits) << value << std::endl;
   }
   void bool_print(const char* label, bool value)
   {  std::cout << std::setw(35) << std::left << label;
      if( value )
         std::cout << " = yes";
      else
         std::cout << " = no";
      std::cout << std::endl;
   }
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
}

int main(int argc, const char* argv[])
{
   bool   ok = true;
   double inf = std::numeric_limits<double>::infinity();
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
      "start_near_solution"
   };
   size_t n_arg = sizeof(arg_name) / sizeof(arg_name[0]);
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
   assert( n_arg == 9 );
   size_t random_seed            = std::atoi( argv[1] );
   size_t number_random          = std::atoi( argv[2] );
   bool   quasi_fixed            = string( argv[3] ) == "yes";
   bool   trace_optimize_fixed   = string( argv[4] ) == "yes";
   bool   ipopt_solve            = string( argv[5] ) == "yes";
   bool   bool_sparsity          = string( argv[6] ) == "yes";
   bool   hold_memory            = string( argv[7] ) == "yes";
   bool   derivative_test        = string( argv[8] ) == "yes";
   bool   start_near_solution    = string( argv[9] ) == "yes";
   //
   // hold memory setting
   CppAD::thread_alloc::hold_memory(hold_memory);
   //
   // print the command line arugments with labels for each value
   for(size_t i = 0; i < n_arg; i++)
      label_print(arg_name[i], argv[1+i]);
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
   // initialize gsl random number generator
   size_t actual_seed = CppAD::mixed::new_gsl_rng(random_seed);
   label_print("actual_seed", actual_seed);
   //
   size_t n_data   = number_random;
   size_t n_random = number_random;
   size_t n_fixed  = 1;
   //
   // explicit constriants
   d_vector fix_constraint_lower(0), fix_constraint_upper(0);
   //
   //
   // difference
   gsl_rng* rng = CppAD::mixed::get_gsl_rng();
   d_vector y(n_data), random_in(n_random), theta_sim(n_fixed);
   theta_sim[0] = 1.0;
   for(size_t i = 0; i < n_data; i++)
   {  y[i]         = double(i + 1) * theta_sim[0];
      y[i]        += gsl_ran_gaussian(rng, sigma_y);
      random_in[i] = 0.0;
   }
   //
   d_vector
      fixed_lower(n_fixed), fixed_in(n_fixed), fixed_upper(n_fixed);
   fixed_lower[0] = 1e-5; fixed_in[0] = 2.0; fixed_upper[0] = inf;
   if( start_near_solution )
      fixed_in[0] = theta_sim[0];
   //
   // object that is derived from cppad_mixed
   CppAD::mixed::d_sparse_rcv A_rcv; // empty matrix
   mixed_derived mixed_object(
      n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv, y
   );
   // initialization
   size_t thread         = CppAD::thread_alloc::thread_num();
   size_t start_bytes    = CppAD::thread_alloc::inuse(thread);
   double start_seconds  = CppAD::elapsed_seconds();
   //
   mixed_object.initialize(fixed_in, random_in);
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
   // optimize the fixed effects using quasi-Newton method
   string random_ipopt_options =
      "Integer print_level               0\n"
      "String  sb                        yes\n"
      "String  derivative_test           none\n"
      "Numeric tol                       1e-7\n"
   ;
   if( ipopt_solve ) random_ipopt_options +=
      "String  evaluation_method         ipopt_solve\n";
   string fixed_ipopt_options =
      "String  sb                        yes\n"
      "Numeric tol                       1e-7\n"
      "Integer max_iter                  15\n"
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
   d_vector random_lower(n_random), random_upper(n_random);
   for(size_t i = 0; i < n_random; i++)
   {  random_lower[i] = -inf;
      random_upper[i] = +inf;
   }
   // optimize fixed effects
   start_seconds = CppAD::elapsed_seconds();
   d_vector fixed_scale = fixed_in;
   CppAD::mixed::fixed_solution solution = mixed_object.optimize_fixed(
      fixed_ipopt_options,
      random_ipopt_options,
      fixed_lower,
      fixed_upper,
      fix_constraint_lower,
      fix_constraint_upper,
      fixed_scale,
      fixed_in,
      random_lower,
      random_upper,
      random_in
   );
   end_seconds = CppAD::elapsed_seconds();
   if( trace_optimize_fixed )
      std::cout << std::endl;
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
      random_lower,
      random_upper,
      random_in
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
   size_t number_fixed_samples = 1000;
   d_vector sample( number_fixed_samples * n_fixed );
   start_seconds = CppAD::elapsed_seconds();
   mixed_object.sample_fixed(
      sample,
      hes_fixed_obj_rcv,
      solution,
      fixed_lower,
      fixed_upper
   );
   end_seconds = CppAD::elapsed_seconds();
   label_print("sample_fixed_seconds", end_seconds - start_seconds);
   //
   end_bytes          = CppAD::thread_alloc::inuse(thread);
   string final_bytes = size_t2string(end_bytes - start_bytes);
   label_print("final_bytes", final_bytes);
   //
   // compute the sample standard deviations
   // and ratio of error divided by sample standard deviation
   d_vector sample_std(n_fixed), estimate_ratio(n_fixed);
   for(size_t j = 0; j < n_fixed; j++)
      sample_std[j] = 0.0;
   for(size_t i = 0; i < number_fixed_samples; i++)
   {  for(size_t j = 0; j < n_fixed; j++)
      {  double diff    = sample[i * n_fixed + j] - theta_out[j];
         sample_std[j] += diff * diff;
      }
   }
   for(size_t j = 0; j < n_fixed; j++)
   {  sample_std[j] = std::sqrt(sample_std[j]/double(number_fixed_samples));
      // check if results results are reasonable
      estimate_ratio[j] = ( theta_out[j] - theta_sim[j] ) / sample_std[j];
      ok  &= std::fabs(estimate_ratio[j]) < 10.0;
   }
   label_print("theta_0_estimate", theta_out[0] );
   label_print("theta_0_std",      sample_std[0] );
   label_print("theta_0_ratio",    estimate_ratio[0] );
   //
   // Check that only the lower limit on theta[1] is active.
   ok &= solution.fixed_lag.size() == n_fixed;
   ok &= solution.fixed_lag[0] == 0.0;
   ok &= solution.fix_con_lag.size() == 0;
   ok &= solution.ran_con_lag.size() == 0;
   //
   if( ok )
      label_print("ar1_xam_ok", "yes");
   else
      label_print("ar1_xam_ok", "no");
   //
   CppAD::thread_alloc::hold_memory(false);
   CppAD::mixed::free_gsl_rng();
   if( ok )
      return 0;
   return 1;
}
// END C++
