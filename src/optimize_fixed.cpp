// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
/*
{xrst_begin optimize_fixed}
{xrst_spell
  iter
  nlp
  struct
}

Optimize Fixed Effects
######################

Syntax
******
*solution*  =

| *mixed_object* . ``optimize_fixed`` (
| |tab| *fixed_ipopt_options* ,
| |tab| *random_ipopt_options* ,
| |tab| *fixed_lower* ,
| |tab| *fixed_upper* ,
| |tab| *fix_constraint_lower* ,
| |tab| *fix_constraint_upper* ,
| |tab| *fixed_scale* ,
| |tab| *fixed_in* ,
| |tab| *random_lower* ,
| |tab| *random_upper* ,
| |tab| *random_in* ,
| |tab| *warm_start*
| )

Purpose
*******
This routine maximizes the fixed effects objective
:ref:`L(theta)<theory@Objective@Fixed Effects Objective, L(theta)>` .

inf
***
The value ``inf`` below refers to

   ``std::numeric_limits<double>::infinity`` ()

mixed_object
************
We use :ref:`derived_ctor@mixed_object`
to denote an object of a class that is
derived from the ``cppad_mixed`` base class.

fixed_ipopt_options
*******************
This argument has prototype

   ``const std::string&`` *fixed_ipopt_options*

and is the :ref:`ipopt_options-name` for optimizing the fixed effects
with the following qualifications:

derivative_test
===============
If *derivative_test* is ``none`` ,
no derivative testing is done.
If it is ``first-order`` ,
only first order derivatives are tested.
If it is ``second-order`` (``only-second-order`` )
second derivatives (and first derivatives) are tested.
In these two cases
:ref:`derived_ctor@quasi_fixed` must be ``false`` .
If *derivative_test* is ``adaptive`` ,
a special ``cppad_mixed`` adaptive step size method is
used to test first order derivatives.
If it is ``trace-adaptive`` ,
the adaptive step size results are traced on standard output.

hessian_approximation
=====================
If *quasi_fixed* is true (false),
*hessian_approximation* will be set to
``limit-memory`` (``exact`` ).
If it is also set in *fixed_ipopt_options* , it must have this value.

max_iter
========
If *max_iter* == ``-1`` in *fixed_ipopt_options* ,
*solution* . ``fixed_opt`` == *fixed_in* .
In addition, Ipopt is run with *max_iter*  = 0 and the return status
``Ipopt::Maximum_Iterations_Exceeded`` is consider normal; i.e.,
does not generate a warning or error message.
Furthermore, the fixed effects optimization will return immediately
(not try to backup and recover) if an error occur during evaluation of
the fixed effects objective, constraints, or their derivatives.
If *max_iter*  == 0 in the options,
it may be that *solution* . ``fixed_opt`` != *fixed_in*
(Ipopt moves non-equality constraints to the interior of the constraint).
(this is the only difference between ``-1`` and ``0`` ).

accept_after_max_steps
======================
The default value for *accept_after_max_steps* is ``-1``
(no limit).
This is the maximum number of backtracking steps to take before
accepting a line search point; see
:ref:`ipopt_trace@ls` is the ipopt tracing documentation.

nlp_scaling_method
==================
When optimizing the fixed effects,
the objective and the constraint functions are automatically scaled
by ``cppad_mixed`` .
It is an error for the user to specify this option in
*fixed_ipopt_options* .

random_ipopt_options
********************
This argument has prototype

   ``const std::string&`` *random_ipopt_options*

and is the :ref:`ipopt_options-name` for optimizing the random effects.

fixed_lower
***********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_lower*

It has size :ref:`derived_ctor@n_fixed` and
specifies the lower limits for the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` .
Note that minus infinity is used for no lower limit.

fixed_upper
***********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_upper*

It has size *n_fixed* and specifies the upper limits for the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>` .
Note that plus infinity is used for no upper limit.

fix_constraint_lower
********************
This argument has prototype

   ``const CppAD::vector<double>&`` *fix_constraint_lower*

it has size *n_fixed* and specifies the lower limits for the
:ref:`fixed constraints<fix_constraint-name>` .
Note that minus infinity is used for no lower limit.

fix_constraint_upper
********************
This argument has prototype

   ``const CppAD::vector<double>&`` *fix_constraint_upper*

it specifies the upper limits for the
:ref:`fixed constraints<fix_constraint-name>` .
Note that plus infinity is used for no upper limit.

fixed_scale
***********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_scale*

The fixed effect objective and constraint functions are multiplied by a
scale factor so that their derivatives are near one at *fixed_scale* .
This makes the Ipopt tolerance be relative to the derivatives at
*fixed_scale* .
It must hold for each *j* that

   *fixed_lower* [ *j* ] <= *fixed_scale* [ *j* ] <= *fixed_upper* [ *j* ]

Partial derivatives with respect to components for which

   *fixed_lower* [ *j* ] == *fixed_upper* [ *j* ]

are not included in this scaling.
Note that you can continue an optimization with the same scaling
by setting

   *fixed_in* = *solution* . ``fixed_opt``

and the re-running the optimization.
Also note that scaling the fixed effects is not done by ``cppad_mixed``
and should be done by the users program when it is useful.

fixed_in
********
This argument has prototype

   ``const CppAD::vector<double>&`` *fixed_in*

It specifies the initial value for the
:ref:`fixed effects<problem@Notation@Fixed Effects, theta>`
vector :math:`\theta` during the optimization process.
It must hold for each *j* that

   *fixed_lower* [ *j* ] <= *fixed_in* [ *j* ] <= *fixed_upper* [ *j* ]

random_lower
************
This argument has prototype

   ``const CppAD::vector<double>&`` *random_lower*

It must have size equal to
:ref:`derived_ctor@n_random` and
specifies the lower limits for the optimization of the
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
This may be useful to keep the random effects
out of regions of numerical instability.
On the other hand, the calculation of the
:ref:`theory@Derivative of Optimal Random Effects`
:math:`\hat{u}_\theta ( \theta )` will not be correct when these constraints
are active (and this could have adverse effects on the optimization).
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
:ref:`random effects<problem@Notation@Random Effects, u>`
vector :math:`u`.
It must hold that

   *random_lower* [ *i* ] <= *random_in* [ *i* ] <= *random_upper* [ *i* ]

for each valid index *i* .

warm_start
**********
This argument is optional and has prototype

   ``const warm_start_struct&`` *warm_start*

It is an error for the user to specify *warm_start_init_point* in
*fixed_ipopt_options* .

No Warm Start
=============
If the size of *warm_start* . ``x_info`` is zero,
there is no warm start information.
This is the same as when the argument is not present.
In this case, the ipopt *warm_start_init_point* option will be set to
``no`` .

Warm Start
==========
If the size of *warm_start* . ``x_info`` is non-zero,
*warm_start* must is equal the
:ref:`fixed_solution@warm_start` field in a
fixed effects solution returned by a previous call to ``optimized_fixed`` .
This can be used to continue a fit when the maximum number of iterations
is reached or when the tolerance for the fixed or random effects is changed.

#. The ipopt *warm_start_init_point* options will be set to ``yes`` .
#. The ipopt *mu_strategy* options will be set to ``monotone`` .
#. The *fixed_scale* and *fixed_in* arguments are not used
   during a warm start optimization.
#. Any derivative test specified in *fixed_ipopt_options*
   will not be passed onto Ipopt; i.e.,
   no derivative testing is done during a warm start.

Example
=======
see :ref:`warm_start.cpp-name`

solution
********
The return value has prototype

   ``CppAD::mixed::fixed_solution`` *solution*

It is the solution (obtained by optimization) of the
fixed effects vector and its Lagrange multipliers; see
:ref:`fixed_solution-name` .

Laplace Approximation
*********************
The :ref:`theory-name` for the
Laplace approximation optimization only includes the case where
the :ref:`random likelihood<ran_likelihood-name>` is smooth.

{xrst_comment ipopt_options is also used by optimize_random}
{xrst_toc_hidden
   example/user/optimize_fixed.cpp
   xrst/ipopt_options.xrst
   xrst/ipopt_trace.xrst
}
Example
*******
The file :ref:`optimize_fixed.cpp-name` contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

ipopt_fixed
***********
The  class ``ipopt_fixed`` is used by ``optimize_fixed``
to optimize the fixed effects.
It's specifications are not part of the :ref:`cppad_mixed-name` public interface.

{xrst_end optimize_fixed}
------------------------------------------------------------------------------
*/
# include <coin-or/IpIpoptApplication.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/ipopt_app_status.hpp>

CppAD::mixed::fixed_solution cppad_mixed::try_optimize_fixed(
   const std::string&     fixed_ipopt_options           ,
   const std::string&     random_ipopt_options          ,
   const d_vector&        fixed_lower                   ,
   const d_vector&        fixed_upper                   ,
   const d_vector&        fix_constraint_lower          ,
   const d_vector&        fix_constraint_upper          ,
   const d_vector&        fixed_scale                   ,
   const d_vector&        fixed_in                      ,
   const d_vector&        random_lower                  ,
   const d_vector&        random_upper                  ,
   const d_vector&        random_in                     ,
   const CppAD::mixed::warm_start_struct& warm_start    )
{  bool ok = true;
   using Ipopt::SmartPtr;
   assert( warm_start.x_info.size() != 0 || warm_start.g_info.size() == 0 );
   //
   // fixed_(lower, upper, in)
   assert( fixed_lower.size() == n_fixed_ );
   assert( fixed_upper.size() == n_fixed_ );
   assert( fixed_in.size()    == n_fixed_ );
   //
   // fix_constraint(lower and upper)
   assert( fix_constraint_lower.size() == fix_con_fun_.Range() );
   assert( fix_constraint_upper.size() == fix_con_fun_.Range() );
   //
   // random_(lower, upper, in)
   assert( random_lower.size() == n_random_ );
   assert( random_upper.size() == n_random_ );
   assert( random_in.size() == n_random_ );
   //
   // make sure initialize has been called
   if( ! initialize_done_ )
   {  std::string error_message =
      "cppad_mixed::initialize was not called before optimize_fixed";
      fatal_error(error_message);
   }
   //
   // optimal random effect for initial fixed effects
   d_vector random_opt;
   if( n_random_ > 0 ) random_opt = try_optimize_random(
      random_ipopt_options,
      fixed_in,
      random_lower,
      random_upper,
      random_in
   );
   if( ! init_laplace_obj_done_ && ! quasi_fixed_ && n_random_ > 0 )
   {  init_laplace_obj(fixed_in, random_opt);
      assert( init_laplace_obj_done_ );
      assert( init_laplace_obj_fun_done_ );
      assert( init_laplace_obj_hes_done_ );
   }
# ifndef NDEBUG
   for(size_t j = 0; j < n_fixed_; j++)
   {  assert( fixed_lower[j] <= fixed_in[j] );
      assert( fixed_in[j]    <= fixed_upper[j] );
   }
# endif
   // create a reference to this object
   cppad_mixed& mixed_object(*this);
   //
   // Create an instance of an IpoptApplication
   SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
   //
   app->Options()->SetStringValue("nlp_scaling_method", "none");
   //
   // hessian_approximation
   if( quasi_fixed_ )
   {  app->Options()->SetStringValue(
         "hessian_approximation", "limited-memory");
   }
   else
   {  app->Options()->SetStringValue(
         "hessian_approximation", "exact");
   }
   // warm_start_init_point
   if( warm_start.x_info.size() == 0 )
   {  app->Options()->SetStringValue( "warm_start_init_point", "no");
   }
   else
   {  app->Options()->SetStringValue( "warm_start_init_point", "yes");
      app->Options()->SetStringValue( "mu_strategy", "monotone");
      app->Options()->SetNumericValue("mu_init", warm_start.mu);
   }
   //
   // accept_after_max_steps
   app->Options()->SetIntegerValue(
      "accept_after_max_steps", -1
   );
   // get the initial value for maximum number of iterations
   std::string tag    = "max_iter";
   std::string prefix = "";
   Ipopt::Index fixed_max_iter;
   app->Options()->GetIntegerValue(tag, fixed_max_iter, prefix);
   assert( fixed_max_iter >= 0 );
   //
   // default for adaptive derivative test
   std::string adaptive_check = "none";
   // ----------------------------------------------------------------------
   // Set options for optimization of the fixed effects
   const std::string& options = fixed_ipopt_options;
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
      // switch on option type
      if ( tok_1 == "String" )
      {  if( tok_2 == "derivative_test" )
         {  if( warm_start.x_info.size() != 0 )
            {  // no derivative testing during a warm start
               adaptive_check = "none";
               tok_3          = "none";
            }
            if( tok_3 == "adaptive" || tok_3 == "trace-adaptive" )
            {  adaptive_check = tok_3;
               tok_3 = "none";
            }
         }
         if( tok_2 == "nlp_scaling_method" )
         {  std::string msg = "optimize_fixed: fixed_ipopt_options: ";
            msg += " nlp_scaling_method cannot be specified";
            fatal_error(msg);
         }
         if( tok_2 == "warm_start_init_point" )
         {  std::string msg = "optimize_fixed: fixed_ipopt_options: ";
            msg += " warm_start_init_point cannot be specified";
            fatal_error(msg);
         }
         app->Options()->SetStringValue(tok_2.c_str(), tok_3.c_str());
         if( quasi_fixed_ )
         {  ok = true;
            if( tok_2 == "hessian_approximation" )
               ok &= tok_3 == "limited-memory";
            if( tok_2 == "derivative_test" )
               ok &= tok_3 == "none" || tok_3 == "first-order";
            if( ! ok )
            {  std::string msg = "cppad_mixed: constructed with";
               msg += " quasi_fixed true so cannot have ";
               msg += tok_2 + " equal to " + tok_3;
               fatal_error(msg);
            }
         }
         else
         {  ok = true;
            if( tok_2 == "hessian_approximation" )
               ok &= tok_3 == "exact";
            if( ! ok )
            {  std::string msg = "cppad_mixed: constructed with";
               msg += " quasi_fixed false so cannot have ";
               msg += tok_2 + " equal to " + tok_3;
               fatal_error(msg);
            }
         }
      }
      else if ( tok_1 == "Numeric" )
      {  Ipopt::Number value = std::atof( tok_3.c_str() );
         app->Options()->SetNumericValue(tok_2.c_str(), value);
      }
      else if ( tok_1 == "Integer" )
      {  Ipopt::Index value = std::atoi( tok_3.c_str() );
         if( tok_2 == "max_iter" )
         {  fixed_max_iter = value;
            if( value == -1 )
               value = 0;
         }
         app->Options()->SetIntegerValue(tok_2.c_str(), value);
      }
      else assert(false);
      //
      // next line
      begin_1 = end_3;
      while( options[begin_1] == ' ' || options[begin_1] == '\n' )
         begin_1++;
   }
   // ----------------------------------------------------------------------
   // when fixed_max_iter is -1 we want to know what is failing at the current
   // values of the fixed effects, so abort right away instead of backing up
   bool abort_on_eval_error = fixed_max_iter == -1;
   // ----------------------------------------------------------------------
   // get the tolerance setting for the fixed effects optimization
   tag    = "tol";
   prefix = "";
   double fixed_tolerance;
   app->Options()->GetNumericValue(tag, fixed_tolerance, prefix);
   //
   // object used to evalutate objective, constraints, and their derivatives
   // (note that  one does not need to delete an ipopt smart pointer)
   SmartPtr<CppAD::mixed::ipopt_fixed> fixed_nlp =
   new CppAD::mixed::ipopt_fixed(
      warm_start,
      abort_on_eval_error,
      random_ipopt_options,
      fixed_tolerance,
      fixed_lower,
      fixed_upper,
      fix_constraint_lower,
      fix_constraint_upper,
      fixed_scale,
      fixed_in,
      random_lower,
      random_upper,
      random_opt,
      mixed_object
   );
   // error message is empty directly after constructor
   assert( fixed_nlp->get_error_message() == "" );
   assert( ok );

   // check derivative calculation ?
   double relative_tol = std::numeric_limits<double>::infinity();
   bool   trace        = false;
   if( adaptive_check != "none" )
   {  if( adaptive_check == "trace-adaptive" )
         trace = true;
      else
      {  assert( adaptive_check == "adaptive" );
      }
      relative_tol  = 1e-3;
   }
   // if trace is false and relative_tol is infinity, no derivative check
   // is done but the scale factors scale_f_ and scale_g_ are still computed.
   ok = fixed_nlp->adapt_derivative_chk(trace, relative_tol);
   if( fixed_nlp->get_error_message() != "" )
   {  std::string msg = "optimize_fixed: ";
      msg            += fixed_nlp->get_error_message();
      fatal_error(msg);
   }
   if( ! ok )
   {  warning("optimize_fixed: adaptive derivative test failed");
   }
   //
   // Set values used for minus and plus infinity
   app->Options()->SetNumericValue(
      "nlp_lower_bound_inf", fixed_nlp->nlp_lower_bound_inf()
   );
   //
   // variable to hold status values returned by app
   Ipopt::ApplicationReturnStatus status;
   //
   // initialize app
   status = app->Initialize();
   ok    &= status == Ipopt::Solve_Succeeded;
   if( ! ok )
   {  std::string message = "optimize_fixed: ";
      message += CppAD::mixed::ipopt_app_status(status);
      if( status == Ipopt::Maximum_Iterations_Exceeded )
         warning(message);
      else
         fatal_error(message);
   }
   //
   // all_constrained
   bool all_constrained = true;
   for(size_t j = 0; j < n_fixed_; ++j)
      all_constrained &= fixed_lower[j] == fixed_upper[j];
   //
   // fixed_nlp.solution_
   if( all_constrained )
   {  // ipopt crashes when all_constrained is true,
      // so just set the solution to fixed_in
      fixed_nlp -> fixed_eq_constrain( fixed_in );
   }
   else
   { // solve the problem
      status = app->OptimizeTNLP(fixed_nlp);
   }
   //
   // check if fixed_nlp could not evaluate one of its functions
   if(  fixed_nlp->get_error_message() != "" )
   {  // report the error message as a warning
      warning( fixed_nlp->get_error_message() );
# ifndef NDEBUG
      if( status == Ipopt::User_Requested_Stop )
         assert( fixed_nlp->get_error_message() != "" );
# endif
      // ---------------------------------------------------------------
      // perhaps the following will generate a more useful error message.
      a1_vector a1_theta(n_fixed_);
      for(size_t j = 0; j < n_fixed_; j++)
         a1_theta[j] = a1_double( fixed_nlp->get_error_fixed(j) );
      try {
         fix_constraint(a1_theta);
         fix_likelihood(a1_theta);
      }
      catch(const std::exception& e)
      {  std::string error_message = "optimize_fixed: std::exception: ";
         error_message += e.what();
         fatal_error(error_message);
      }
      catch(const CppAD::mixed::exception& e)
      {  std::string error_message = e.message("optimize_fixed");
         warning(error_message);
      }
      // ---------------------------------------------------------------
   }

   //
   if( fixed_max_iter <= 0 )
   {  // special case where we are not trying to solve
      assert( ok );
      ok  = status == Ipopt::Maximum_Iterations_Exceeded;
      ok |= status == Ipopt::Solve_Succeeded;
      if( ! ok )
      {  warning("optimize_fixed: unexpected error during zero iterations");
      }
   }
   else if( ! all_constrained )
   {
      if( status != Ipopt::Solve_Succeeded )
      {  warning("optimize_fixed: ipopt failed to converge");
      }
      if( ! fixed_nlp->finalize_solution_ok_ )
      {  warning(
            "optimize_fixed: double check of Ipopt solution failed\n"
            "bounds, complementarity, or gradient of Lagragian test."
         );
      }
   }
   //
   // return the entire solution including the lagrange multipliers
   CppAD::mixed::fixed_solution solution = fixed_nlp->solution();
   //
   // special case were we return fixed_in
   if( fixed_max_iter == -1 )
      solution.fixed_opt = fixed_in;
   //
   return solution;
}
// ---------------------------------------------------------------------------
CppAD::mixed::fixed_solution cppad_mixed::optimize_fixed(
   const std::string& fixed_ipopt_options     ,
   const std::string& random_ipopt_options    ,
   const d_vector&    fixed_lower             ,
   const d_vector&    fixed_upper             ,
   const d_vector&    fix_constraint_lower    ,
   const d_vector&    fix_constraint_upper    ,
   const d_vector&    fixed_scale             ,
   const d_vector&    fixed_in                ,
   const d_vector&    random_lower            ,
   const d_vector&    random_upper            ,
   const d_vector&    random_in               ,
   const CppAD::mixed::warm_start_struct& warm_start )
{  CppAD::mixed::fixed_solution ret;
   try
   {  ret = try_optimize_fixed(
         fixed_ipopt_options     ,
         random_ipopt_options    ,
         fixed_lower             ,
         fixed_upper             ,
         fix_constraint_lower    ,
         fix_constraint_upper    ,
         fixed_scale             ,
         fixed_in                ,
         random_lower            ,
         random_upper            ,
         random_in               ,
         warm_start
      );
   }
   catch(const std::exception& e)
   {  std::string error_message = "optimize_fixed: std::exception: ";
      error_message += e.what();
      fatal_error(error_message);
      assert(false);
   }
   catch(const CppAD::mixed::exception& e)
   {  std::string error_message = e.message("optimize_fixed");
      fatal_error(error_message);
      assert(false);
   }
   return ret;
}
