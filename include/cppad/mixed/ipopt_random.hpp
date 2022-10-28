// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_IPOPT_RANDOM_HPP
# define CPPAD_MIXED_IPOPT_RANDOM_HPP

/*
$begin ipopt_random$$
$spell
   CppAD
   ran_obj
   cppad
   obj
   Ipopt
   nlp
   inf
   hpp
   std
$$

$section Ipopt NLP Class Used to Optimize Random Effects$$

$head Syntax$$
$code # include <cppad/mixed/ipopt_random.hpp>$$

$head Private$$
This class is an implementation detail and not part of the
CppAD Mixed user API.

$head get_error_message()$$
This function returns a $code std::string$$ value that
corresponds to the most recent error message.
If it is non-empty, there is no error.

$head clear_error_message()$$
This function sets the current error message to the empty string.

$head nlp_lower_bound_inf()$$
This member function returns the $code double$$ value used
for minus infinity as a lower bound.

$head nlp_upper_bound_inf()$$
This member function returns the $code double$$ value used
for plus infinity as an upper bound.

$head solution()$$
This member function returns a
$cref fixed_solution$$ object containing the
optimal solution information.

$head Default Destructor$$
The default destructor is defined by this include file.


$childtable%src/ipopt_random/ctor.cpp
   %src/ipopt_random/get_nlp_info.cpp
   %src/ipopt_random/get_bounds_info.cpp
   %src/ipopt_random/get_starting_point.cpp
   %src/ipopt_random/eval_f.cpp
   %src/ipopt_random/eval_grad_f.cpp
   %src/ipopt_random/eval_g.cpp
   %src/ipopt_random/eval_jac_g.cpp
   %src/ipopt_random/eval_h.cpp
   %src/ipopt_random/finalize_solution.cpp
%$$

$end
-----------------------------------------------------------------------------
*/
# include <coin-or/IpTNLP.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/fixed_solution.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
   //
   // ipopt_random
   class ipopt_random : public Ipopt::TNLP
   {
   private:
      // cppad_mixed types used by this class
      typedef cppad_mixed::d_vector     d_vector;
      typedef CppAD::vector<size_t>     s_vector;
      //
      // Ipopt types used by this class
      typedef Ipopt::Number               Number;
      typedef Ipopt::Index                Index;
      typedef Ipopt::TNLP::IndexStyleEnum IndexStyleEnum;
      // ---------------------------------------------------------------
      // member variables set during constructor
      size_t          n_fixed_;       // number of fixed effects
      size_t          n_random_;      // number of random effects;
      const d_vector& fixed_vec_;     // value of fixed effects
      const d_vector& random_lower_;  // lower limit for random effects
      const d_vector& random_upper_;  // upper limit for random effects
      const d_vector& random_in_;     // initial value for random effects
      const size_t nnz_h_lag_;        // # non-zeros in Hessian of Lagragian
      // effectively const (not changed after constructor)
      double nlp_lower_bound_inf_;    // Ipopt's code for - infinity
      double nlp_upper_bound_inf_;    // Ipopt's code for + infinity
      // values initialized by constructor
      cppad_mixed&    mixed_object_;  // cppad_mixed for this problem
      // objective corresponding to the most recent x
      double objective_current_;
      // ------------------------------------------------------------------
      // If empty, no error has been detected by the ipopt_random class
      // Otherwise, this was the last error detected.
      // Initialize by constructor to have size n_random_.
      std::string error_message_;
      // ------------------------------------------------------------------
      // if error_message_ is non-empty, this is the most recent
      // random effects that resulted in a CppAD::mixed::exception
      d_vector error_random_;
      // ------------------------------------------------------------------
      // final solution returned by optimization
      d_vector random_opt_;
      // ------------------------------------------------------------------
      // Public versions wrap these functions in a try / catch block.
      // If an eval fails, it sets this message and returns false.
      void try_eval_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number&         obj_value
      );
      void try_eval_grad_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number*         grad_f
      );
      void try_eval_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Number*         g
      );
      void try_eval_jac_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Index           nele_jac ,
         Index*          iRow     ,
         Index*          jCol     ,
         Number*         values
      );
      void try_eval_h(
         Index         n              ,
         const Number* x              ,
         bool          new_x          ,
         Number        obj_factor     ,
         Index         m              ,
         const Number* lambda         ,
         bool          new_lambda     ,
         Index         nele_hess      ,
         Index*        iRow           ,
         Index*        jCol           ,
         Number*       values
      );
   public:
      //  get the current ipopt_random error message
      std::string get_error_message(void) const
      {  return error_message_; }
      //
      //  get the random effects corresponding to the error message
      double get_error_random(size_t j) const
      {  return error_random_[j]; }
      //
      //  clear the current ipopt_random error message
      void clear_error_message(void)
      {  error_message_ = ""; }
      //
      // get minus infinity
      double nlp_lower_bound_inf(void) const
      {  return nlp_lower_bound_inf_; }
      //
      // get plus infinity
      double nlp_upper_bound_inf(void) const
      {  return nlp_upper_bound_inf_; }
      //
      // optimal solution
      d_vector random_opt(void) const
      {  return random_opt_; }
      // -----------------------------------------------------------------
      // constructor
      ipopt_random(
         const d_vector& fixed_vec    ,
         const d_vector& random_lower ,
         const d_vector& random_upper ,
         const d_vector& random_in    ,
         cppad_mixed&    mixed_object
      );
      //
      // destructor
      virtual ~ipopt_random(void)
      { }
      // -----------------------------------------------------------------
      virtual bool get_nlp_info(
         Index&          n            ,
         Index&          m            ,
         Index&          nnz_jac_g    ,
         Index&          nnz_h_lag    ,
         IndexStyleEnum& index_style
      );
      virtual bool get_bounds_info(
            Index       n        ,
            Number*     x_l      ,
            Number*     x_u      ,
            Index       m        ,
            Number*     g_l      ,
            Number*     g_u
      );
      virtual bool get_starting_point(
         Index           n            ,
         bool            init_x       ,
         Number*         x            ,
         bool            init_z       ,
         Number*         z_L          ,
         Number*         z_U          ,
         Index           m            ,
         bool            init_lambda  ,
         Number*         lambda
      );
      virtual bool eval_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number&         obj_value
      );
      virtual bool eval_grad_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number*         grad_f
      );
      virtual bool eval_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Number*         g
      );
      virtual bool eval_jac_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Index           nele_jac ,
         Index*          iRow     ,
         Index*          jCol     ,
         Number*         values
      );
      virtual bool eval_h(
         Index         n              ,
         const Number* x              ,
         bool          new_x          ,
         Number        obj_factor     ,
         Index         m              ,
         const Number* lambda         ,
         bool          new_lambda     ,
         Index         nele_hess      ,
         Index*        iRow           ,
         Index*        jCol           ,
         Number*       values
      );
      virtual void finalize_solution(
         Ipopt::SolverReturn               status    ,
         Index                             n         ,
         const Number*                     x         ,
         const Number*                     z_L       ,
         const Number*                     z_U       ,
         Index                             m         ,
         const Number*                     g         ,
         const Number*                     lambda    ,
         Number                            obj_value ,
         const Ipopt::IpoptData*           ip_data   ,
         Ipopt::IpoptCalculatedQuantities* ip_cq
      );
   };
} } // END_CPPAD_MIXED_NAMESPACE

# endif
