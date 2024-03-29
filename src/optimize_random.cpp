// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/ipopt/solve.hpp>
# include <coin-or/IpIpoptApplication.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/ipopt_random.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_app_status.hpp>
/*
{xrst_begin optimize_random}

Optimize Random Effects
#######################

Syntax
******
*random_out*  =

| *mixed_object* . ``optimize_random`` (
| |tab| *options* , *fixed_vec* , *random_lower* , *random_upper* , *random_in*
| )

Purpose
*******
This routine maximizes the
:ref:`random likelihood<ran_likelihood-name>`
corresponding to the object *mixed_object* .

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

options
*******
This argument  has the prototype

   ``const std::string&`` *options*

This is the :ref:`ipopt_options-name` for optimizing the random effects
with the following qualifications:

evaluation_method
=================
There is an additional :ref:`ipopt_options@String` option with
*name* = ``evaluation_method``
and *value* is either ``ipopt_random``  or ``ipopt_solve`` .
The ``ipopt_random`` choice uses
``CppAD::mixed::ipopt_random`` for optimizing random effects.
This special purpose class
is expected to eventually be the faster choice.
This is the default choice; i.e., ``ipopt_random`` will be used
if this option is not present.
The ``ipopt_solve`` choice uses
``CppAD::ipopt::solve`` for optimizing random effects.
Currently this is sometimes faster and so this choice is still included
(but may be removed in the future).

fixed_vec
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_vec*

It specifies the value of the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta`.

random_lower
************
This argument has prototype

   ``const CppAD::vector<double>&`` *random_lower*

It must have size equal to
:ref:`derived_ctor@n_random` and
specifies the lower limits for the optimization of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
The value minus infinity can be used to specify no lower limit.

random_upper
************
This argument has prototype

   ``const CppAD::vector<double>&`` *random_upper*

It must have size equal to
:ref:`derived_ctor@n_random` and
specifies the upper limits for the optimization of the random effect.
The value plus infinity can be used to specify no lower limit.

random_in
*********
This argument has prototype

   ``const CppAD::vector<double>&`` *random_in*

It must have size equal to
:ref:`derived_ctor@n_random` and
specifies the initial value used for the optimization of the
:ref:`random effects<problem@Notation@Random Effects, u>` vector :math:`u`.
It must hold that

   *random_lower* [ *i* ] <= *random_in* [ *i* ] <= *random_upper* [ *i* ]

for each valid index *i* .

random_out
**********
The return value  has prototype

   ``CppAD::vector<double>`` *random_out*

It is the final value (obtained by optimization) of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
{xrst_toc_hidden
   example/user/optimize_random.cpp
}
Example
*******
The file :ref:`optimize_random.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

{xrst_end optimize_random}
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
   // ------------------------------------------------------------------------
   // class used by use_ipopt_solve
   // (not in empty namespace so can be friend of cppad_mixed class)
   class optimize_random_ipopt {
   public:
      // same as CppAD::mixed::d_vector
      typedef CppAD::vector<double>              Dvector;

      // same as CppAD::mixed::a1_vector
      typedef CppAD::vector< CppAD::AD<double> > ADvector;
   private:
      const size_t     n_fixed_;
      const size_t     n_random_;
      ADvector         fixed_vec_;
      cppad_mixed&     mixed_object_;
   public:
      // constructor
      optimize_random_ipopt(
         size_t           n_random        ,
         const Dvector&   fixed_vec       ,
         cppad_mixed&     mixed_object
      ) :
      n_fixed_      ( fixed_vec.size() ) ,
      n_random_     ( n_random )         ,
      mixed_object_( mixed_object )
      {  fixed_vec_.resize( n_fixed_ );
         for(size_t i = 0; i < n_fixed_; i++)
            fixed_vec_[i] = fixed_vec[i];
      }
      //
      // evaluate objective and constraints
      void operator()(ADvector& fg, const ADvector& x)
      {  assert( fg.size() == 1 );
         //
         // extract the random effects from x
         ADvector random_vec(n_random_);
         for(size_t j = 0; j < n_random_; j++)
            random_vec[j] = x[j];
         //
         // compute log-density vector
         ADvector both_vec(n_fixed_ + n_random_);
         mixed_object_.pack(fixed_vec_, random_vec, both_vec);
         ADvector vec = mixed_object_.ran_like_a1fun_.Forward(0, both_vec);
         assert( vec.size() == 1 );
         //
         // return negative log-likelihood
         fg[0] = vec[0];
      }
   };
} } // END_CPPAD_MIXED_NAMESPACE
namespace { // BEGIN_EMPTY_NAMESPACE
   // ---------------------------------------------------------------------
   // set_options
   std::string set_options(
      Ipopt::SmartPtr<Ipopt::IpoptApplication>& app               ,
      const std::string&                       options            ,
      std::string&                             evaluation_method  )
   {  Ipopt::SmartPtr<Ipopt::IpoptApplication> null(NULL);
      //
      // initilaize return value
      std::string ret;
      //
      // Set options for optimization of the random effects
      size_t begin_1, end_1, begin_2, end_2, begin_3, end_3;
      begin_1     = 0;
      while( options[begin_1] == ' ')
         begin_1++;
      while( begin_1 < options.size() )
      {  // split this line into tokens
         end_1   = options.find_first_of(" \n", begin_1);
         begin_2 = end_1;
         while( options[begin_2] == ' ')
            begin_2++;
         end_2   = options.find_first_of(" \n", begin_2);
         begin_3 = end_2;
         while( options[begin_3] == ' ')
            begin_3++;
         end_3   = options.find_first_of(" \n", begin_3);

         // check for three non-empty tokens
         assert( end_3 != std::string::npos );
         assert( begin_1 < end_1 && end_1 <= begin_2 );
         assert( begin_2 < end_2 && end_2 <= begin_3 );
         assert( begin_3 < end_3 );

         // get the three tokens
         std::string tok_1 = options.substr(begin_1, end_1 - begin_1);
         std::string tok_2 = options.substr(begin_2, end_2 - begin_2);
         std::string tok_3 = options.substr(begin_3, end_3 - begin_3);
         //
         // switch on option type
         bool include = false; // include in app settings
         if ( tok_1 == "String" )
         {  if( tok_2 == "evaluation_method" )
               evaluation_method = tok_3;
            else
            {  if( app != null ) app->Options()->
                  SetStringValue(tok_2.c_str(), tok_3.c_str());
               include = true;
            }
         }
         else if ( tok_1 == "Numeric" )
         {  Ipopt::Number value = std::atof( tok_3.c_str() );
            if( app != null ) app->Options()->
               SetNumericValue(tok_2.c_str(), value);
            include = true;
         }
         else if ( tok_1 == "Integer" )
         {  Ipopt::Index value = std::atoi( tok_3.c_str() );
            if( app != null ) app->Options()->
               SetIntegerValue(tok_2.c_str(), value);
            include = true;
         }
         else assert(false);
         //
         if( include && app == null )
            ret += tok_1 + " " + tok_2 + " " + tok_3 + "\n";
         //
         // next line
         begin_1 = end_3;
         while( options[begin_1] == ' ' || options[begin_1] == '\n' )
            begin_1++;
      }
      return ret;
   }
   // ----------------------------------------------------------------------
   // use_ipopt_random
   CppAD::vector<double> use_ipopt_random(
      const std::string&              options         ,
      const CppAD::vector<double>&    fixed_vec       ,
      const CppAD::vector<double>&    random_lower    ,
      const CppAD::vector<double>&    random_upper    ,
      const CppAD::vector<double>&    random_in       ,
      cppad_mixed&                    mixed_object    )
   {
      // Create an instance of an IpoptApplication
      Ipopt::SmartPtr<Ipopt::IpoptApplication> app =
         IpoptApplicationFactory();
      //
      // Set options for optimization of the random effects
      std::string not_used;
      set_options(app, options, not_used);
      //
      // object that is used to evalute objective and its derivatives
      // (note that one does not need to delete an ipopt smart pointer)
      Ipopt::SmartPtr<CppAD::mixed::ipopt_random> random_nlp =
      new CppAD::mixed::ipopt_random(
         fixed_vec, random_lower, random_upper, random_in, mixed_object
      );
      // error message is empty directory after constructor
      assert( random_nlp->get_error_message() == "" );
      //
      // variable to hold status values returned by app
      Ipopt::ApplicationReturnStatus status;
      //
      // initialize app
      status = app->Initialize();
      if(status != Ipopt::Solve_Succeeded )
         mixed_object.fatal_error("optimize_random: initalization failed");
      //
      // solve the problem
      status = app->OptimizeTNLP(random_nlp);
      //
      // Check if random_nlp could not evaluate one of its functions.
      if(  random_nlp->get_error_message() != "" )
      {  // report the error message as a warning
         mixed_object.warning( random_nlp->get_error_message() );
         //
         typedef CppAD::mixed::a1_vector   a1_vector;
         //
         size_t n_fixed  = fixed_vec.size();
         size_t n_random = random_in.size();
         //
         // Perhaps ran_likelihood will generate a more useful error message
         a1_vector a1_theta(n_fixed), a1_u(n_random);
         for(size_t i = 0; i < n_fixed; i++)
            a1_theta[i] = fixed_vec[i];
         for(size_t i = 0; i < n_random; i++)
            a1_u[i] = random_nlp->get_error_random(i);
         try {
            mixed_object.ran_likelihood(a1_theta, a1_u);
         }
         catch(const std::exception& e)
         {  std::string error_message = "optimize_random: std::exception: ";
            error_message += e.what();
            mixed_object.fatal_error(error_message);
         }
         catch(const CppAD::mixed::exception& e)
         {  std::string error_message = e.message("optimize_random");
            mixed_object.warning(error_message);
         }
      }
      //
      if( status != Ipopt::Solve_Succeeded )
      {  std::string message = "optimize_random: ";
         message += CppAD::mixed::ipopt_app_status(status);
         if( status == Ipopt::Maximum_Iterations_Exceeded )
         {  // In this case, we wish to see how far we got.
            mixed_object.warning(message);
         }
         else if( status == Ipopt::Invalid_Number_Detected )
         {  // Ipopt may move to a point where the objective can be
            // evaluated, but the Hessian can't be evaluated.
            mixed_object.warning(message);
         }
         else throw CppAD::mixed::exception("use_ipopt_random", message);
      }
      //
      return random_nlp->random_opt();
   }
   // ------------------------------------------------------------------------
   // use_ipopt_solve
   CppAD::vector<double> use_ipopt_solve(
      const std::string&              options         ,
      const CppAD::vector<double>&    fixed_vec       ,
      const CppAD::vector<double>&    random_lower    ,
      const CppAD::vector<double>&    random_upper    ,
      const CppAD::vector<double>&    random_in       ,
      cppad_mixed&                    mixed_object    )
   {  typedef CppAD::vector<double> d_vector;
      //
      // infinity
      double inf = std::numeric_limits<double>::infinity();
      //
      // number of independent variable is number of random effects
      // plus number of log-density terms that require absolute values
      size_t nx = random_lower.size();
      //
      // set initial x vector
      d_vector xi(nx);
      for(size_t j = 0; j < nx; j++)
         xi[j] = random_in[j];
      //
      // ipopts default value for infinity (use options to change it)
      double ipopt_infinity = 1e19;
      //
      // set lower and upper limits for x
      d_vector xl(nx), xu(nx);
      for(size_t j = 0; j < nx; j++)
      {  if( random_lower[j] == - inf )
            xl[j] = - ipopt_infinity;
         else
            xl[j] = random_lower[j];
         if( random_lower[j] == + inf )
            xu[j] = + ipopt_infinity;
         else
            xu[j] = random_upper[j];
      }
      //
      // construct fg_eval  object
      CppAD::mixed::optimize_random_ipopt fg_eval(
         nx, fixed_vec, mixed_object
      );
      //
      // optimizer options
      assert( options[ options.size() - 1] == '\n' );
      std::string solve_options = options + "Sparse  true  reverse \n";
      //
      // return solution
      CppAD::ipopt::solve_result<d_vector> solution;
      //
      // empty vectors
      d_vector gl(0), gu(0);
      //
      // solve the optimization problem
      CppAD::ipopt::solve(
         solve_options, xi, xl, xu, gl, gu, fg_eval, solution
      );
      //
      return solution.x;
   }
} // END_EMPTY_NAMESPACE
// ----------------------------------------------------------------------------
// optimize_random
CppAD::vector<double> cppad_mixed::try_optimize_random(
   const std::string& options         ,
   const d_vector&    fixed_vec       ,
   const d_vector&    random_lower    ,
   const d_vector&    random_upper    ,
   const d_vector&    random_in       )
{
   // make sure initialize has been called
   if( ! initialize_done_ )
   {  std::string error_message =
      "optimize_random: initialize was not called before optimize_random";
      fatal_error(error_message);
   }
   if( ! init_ran_like_done_ )
   {  std::string error_message =
      "optimize_random: there are no random effects";
      fatal_error(error_message);
   }
   // number of fixed and random effects
   assert( n_fixed_  == fixed_vec.size() );
   assert( n_random_ == random_in.size() );
   assert( n_random_ == random_lower.size() );
   assert( n_random_ == random_upper.size() );
   //
   // bounds on random effects
# ifndef NDEBUG
   for(size_t j = 0; j < n_random_; j++)
   {  assert( random_lower[j] <= random_in[j] );
      assert( random_in[j]    <= random_upper[j] );
   }
# endif
   // ran_likelihood has no absolute value terms
   assert( ran_like_fun_.Range() == 1 );
   //
   // determine which evaluation method we are using
   std::string evaluation_method = "ipopt_random";
   Ipopt::SmartPtr<Ipopt::IpoptApplication> app = NULL;
   std::string include_options = set_options(app, options, evaluation_method);
   //
   // return value
   d_vector random_opt;
   //
   // this mixed_object
   cppad_mixed& mixed_object(*this);
   //
   if( evaluation_method == "ipopt_solve" )
   {  random_opt = use_ipopt_solve(
         include_options,
         fixed_vec,
         random_lower,
         random_upper,
         random_in   ,
         mixed_object
      );
   }
   else if( evaluation_method == "ipopt_random" )
   {  random_opt = use_ipopt_random(
         include_options,
         fixed_vec,
         random_lower,
         random_upper,
         random_in   ,
         mixed_object
      );
   }
   else
   {  std::string message = "optimize_random: evaluation_method = ";
      message            += evaluation_method + "\n";
      message            += "is not ipopt_solve or ipopt_random";
      fatal_error( message );
   }
   return random_opt;
}
// ----------------------------------------------------------------------------
CppAD::vector<double> cppad_mixed::optimize_random(
   const std::string& options         ,
   const d_vector&    fixed_vec       ,
   const d_vector&    random_lower    ,
   const d_vector&    random_upper    ,
   const d_vector&    random_in       )
{  CppAD::vector<double> ret;
   try
   {  ret = try_optimize_random(
         options, fixed_vec, random_lower, random_upper, random_in
      );
   }
   catch(const std::exception& e)
   {  std::string error_message = "optimize_random: std::exception: ";
      error_message += e.what();
      fatal_error(error_message);
      assert(false);
   }
   catch(const CppAD::mixed::exception& e)
   {  std::string error_message = e.message("optimize_random");
      fatal_error(error_message);
      assert(false);
   }
   return ret;
}
