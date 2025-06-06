// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_IPOPT_FIXED_HPP
# define CPPAD_MIXED_IPOPT_FIXED_HPP

/*
{xrst_begin ipopt_fixed dev}
{xrst_spell
  nlp
}

Ipopt NLP Class Used to Optimize Fixed Effects
##############################################

Private
*******
This class is an implementation detail and not part of the
CppAD Mixed user API.

get_error_message()
*******************
This function returns a ``std::string`` value that
corresponds to the most recent error message.
If it is non-empty, there is no error.

clear_error_message()
*********************
This function sets the current error message to the empty string.

nlp_lower_bound_inf()
*********************
This member function returns the ``double`` value used
for minus infinity as a lower bound.

nlp_upper_bound_inf()
*********************
This member function returns the ``double`` value used
for plus infinity as an upper bound.

solution()
**********
This member function returns a
:ref:`fixed_solution-name` object containing the
optimal solution information.

Default Destructor
******************
The default destructor is defined by this include file.

Contents
********
{xrst_toc_table
   src/ipopt_fixed/ctor.cpp
   src/ipopt_fixed/get_nlp_info.cpp
   src/ipopt_fixed/get_bounds_info.cpp
   src/ipopt_fixed/get_starting_point.cpp
   src/ipopt_fixed/eval_f.cpp
   src/ipopt_fixed/eval_grad_f.cpp
   src/ipopt_fixed/eval_g.cpp
   src/ipopt_fixed/eval_jac_g.cpp
   src/ipopt_fixed/eval_h.cpp
   src/ipopt_fixed/finalize_solution.cpp
   src/ipopt_fixed/fixed_eq_constrain.cpp
   src/ipopt_fixed/intermediate_callback.cpp
   src/ipopt_fixed/adapt_derivative_chk.cpp
   src/ipopt_fixed/one_dim_function.cpp
   src/ipopt_fixed/new_random.cpp
   src/ipopt_fixed/set_scaling.cpp
   example/ipopt_xam.xrst
}

{xrst_end ipopt_fixed}
-----------------------------------------------------------------------------
*/
# include <coin-or/IpTNLP.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/fixed_solution.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
   //
   // ipopt_fixed
   class ipopt_fixed : public Ipopt::TNLP
   {
   private:
      // Ipopt types used by this class
      typedef Ipopt::Number               Number;
      typedef Ipopt::Index                Index;
      typedef Ipopt::TNLP::IndexStyleEnum IndexStyleEnum;
      typedef Ipopt::AlgorithmMode               AlgorithmMode;
      typedef Ipopt::IpoptData                   IpoptData;
      typedef Ipopt::IpoptCalculatedQuantities   IpoptCalculatedQuantities;
      // ---------------------------------------------------------------
      // member variables set during constructor
      //
      const warm_start_struct  warm_start_;
      const bool               abort_on_eval_error_;
      const std::string&       random_ipopt_options_;
      const double             fixed_tolerance_; // ipopt relative tolerance
      //
      const size_t n_fixed_;            // number of fixed effects
      const size_t n_random_;           // number of random effects
      const size_t n_fix_con_;          // number of fixed constraints
      const size_t n_ran_con_;          // number of random constraints
      //
      const d_vector& fixed_lower_;         // fixed effects lower limits
      const d_vector& fixed_upper_;         // fixed effects upper limit
      const d_vector& fix_constraint_lower_;// constraint lower limits
      const d_vector& fix_constraint_upper_;// constraint upper limit
      const d_vector& fixed_scale_;         // fixed effects scalling value
      const d_vector& fixed_in_;            // fixed effects initial value
      //
      const d_vector& random_lower_;    // lower limit for random effects
      const d_vector& random_upper_;    // upper limit for random effects
      const d_vector& random_in_;       // random effects initial value
      //
      cppad_mixed&   mixed_object_;     // cppad_mixed for this problem
      // ---------------------------------------------------------------
      // set during constructor, otherwise const
      double nlp_lower_bound_inf_; // Ipopt's code for - infinity
      double nlp_upper_bound_inf_; // Ipopt's code for + infinity
      //
      size_t fix_likelihood_nabs_; // number of absolute values in prior
      size_t nnz_jac_g_;   // number non-zeros in Jacobian of constraints
      size_t nnz_h_lag_;   // number non-zeros in Hessian of Lagragian
      //
      // maps jac_g sparsity index to row and column index in g'(x).
      s_vector jac_g_row_;
      s_vector jac_g_col_;
      //
      // random constraint jacobian
      CppAD::mixed::d_sparse_rcv ran_con_jac_rcv_;
      //
      // hessian of Laplace objective and constraints
      // (only defined when n_random_ > 0)
      CppAD::mixed::sparse_mat_info laplace_obj_hes_info_;
      //
      s_vector lag_hes_row_;   // row indices for Hessian of Lagrangian
      s_vector lag_hes_col_;   // column indices for Hessian of Lagrangian
      s_vector laplace_obj_hes_2_lag_; // laplace_obj_hes_row_ -> lag_hes_row_
      s_vector fix_like_hes_2_lag_; // fix_like_hes_row_ -> lag_hes_row_
      s_vector fix_con_hes_2_lag_; // fix_con_hes_row -> lag_hes_row
      // ---------------------------------------------------------------
      // This variable is set false by constructor, true at beginning of
      // adapt_derivative_chk
      bool adaptive_called_;
      // The rest of these variables only set by adapt_derivative_chk
      //
      // scale factor (multiplies) components of f(x)
      double scale_f_;
      //
      // scale factor for components of g(x)
      d_vector scale_g_;
      //
      // scale factor for components of x
      d_vector scale_x_;
      // ---------------------------------------------------------------
      // temporaries (size set by constructor only)
      //
      // this vector has size n = n_fixed_ + fix_likelihood_nabs_
      d_vector        x_tmp_;
      //
      d_vector        fixed_tmp_;         // size n_fixed_
      d_vector        c_vec_tmp_;         // size n_fix_con_
      d_vector        A_uhat_tmp_;        // size n_ran_con_
      d_vector        H_beta_tmp_;        // size n_fixed_
      d_vector        w_fix_con_tmp_;     // size n_fix_con_
      d_vector        w_laplace_obj_tmp_; // size n_ran_con_ + 1
      //
      // this vector has size fix_likelihood_nabs_ + 1
      d_vector        w_fix_likelihood_tmp_;
      //
      // if mixed_object_.fix_like_eval returns a vector of size 0, this
      // vector has size 0, otherwise it has size fix_likelihood_nabs_ + 1
      d_vector        fix_likelihood_vec_tmp_;
      // ---------------------------------------------------------------
      // solution to the fixed effects optimization problem
      CppAD::mixed::fixed_solution solution_;
      // ---------------------------------------------------------------
      // Optimal random effects cooressponding to current fixed effects.
      // Set by any eval routine when new_x is true, i.e. new fixed effects.
      d_vector random_cur_;
      // ---------------------------------------------------------------
      // routine used the set random_cur_ (and update ldlt factor)
      void new_random(const d_vector& fixed_vec);
      // ------------------------------------------------------------------
      // If empty, no error has been detected by the ipopt_fixed class
      // Otherwise, this was the last error detected.
      std::string error_message_;
      // ------------------------------------------------------------------
      // if error_message_ is non-empty, this is the most recent
      // fixed effects that threw a CppAD::mixed::exception
      d_vector error_fixed_;
      // ------------------------------------------------------------------
      // Public versions wrap these functions in a try / catch block.
      // If an eval fails, it sets this message and returns false.
      void try_eval_f(
         Index           n        ,
         const d_vector& x        ,
         bool            new_x    ,
         Number&         obj_value
      );
      void try_eval_grad_f(
         Index           n        ,
         const d_vector& x        ,
         bool            new_x    ,
         Number*         grad_f
      );
      void try_eval_g(
         Index           n        ,
         const d_vector& x        ,
         bool            new_x    ,
         Index           m        ,
         Number*         g
      );
      void try_eval_jac_g(
         Index           n        ,
         const d_vector& x        ,
         bool            new_x    ,
         Index           m        ,
         Index           nele_jac ,
         Index*          iRow     ,
         Index*          jCol     ,
         Number*         values
      );
      void try_eval_h(
         Index           n              ,
         const d_vector& x            ,
         bool            new_x          ,
         Number          obj_factor     ,
         Index           m              ,
         const Number*   lambda         ,
         bool            new_lambda     ,
         Index           nele_hess      ,
         Index*          iRow           ,
         Index*          jCol           ,
         Number*         values
      );
      // -------------------------------------------------------------------
      // Used by adapt_derivaive_chk to set scale_f_ and scale_g_
      bool set_scaling(
         const d_vector& x_scale ,
         const d_vector& x_lower ,
         const d_vector& x_upper ,
         const d_vector& grad_f  ,
         const d_vector& jac_g
      );
      // -------------------------------------------------------------------
      // only used by finalize solution
      bool check_in_limits(
         double lower, double x, double upper, double tol
      );
      // -------------------------------------------------------------------
      // Used by adapt_derivative_chk member function to pass
      // information to one_dim_function member function
      //
      // which function is being evaluted along j-th component of x
      enum {
         eval_f_enum,       // f(x)
         eval_g_enum,       // g(x)
         eval_grad_L_enum , // L'(x)
      } one_dim_function_eval_;
      //
      // argument index along which the one dimensional funciton is defined
      size_t      one_dim_function_j_;
      //
      // vector containing all argments
      d_vector    one_dim_function_x_;
      //
      // function factor and Lagrange multipliers used for testing
      double      one_dim_function_obj_factor_;
      d_vector    one_dim_function_lambda_;
   public:
      bool one_dim_function(double x_in, d_vector& fun_out);
      // --------------------------------------------------------------------
      //
      //  get and clear the current ipopt_fixed error message
      std::string get_error_message(void) const
      {  return error_message_; }
      //
      // get the fixed effects corresponding to theh error message
      double get_error_fixed(size_t j) const
      {  return error_fixed_[j]; }
      //
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
      CppAD::mixed::fixed_solution solution(void) const
      {  return solution_; }
      //
      // default destructor
      virtual ~ipopt_fixed(void)
      { }
      // -----------------------------------------------------------------
      //
      // did finalize_solution agree that the solution had converged
      bool finalize_solution_ok_;
      //
      // constructor
      ipopt_fixed(
         const warm_start_struct&  warm_start           ,
         const bool&               abort_on_eval_error  ,
         const std::string&        random_ipopt_options ,
         const double&             fixed_tolerance      ,
         const d_vector&           fixed_lower          ,
         const d_vector&           fixed_upper          ,
         const d_vector&           fix_constraint_lower ,
         const d_vector&           fix_constraint_upper ,
         const d_vector&           fixed_scale          ,
         const d_vector&           fixed_in             ,
         const d_vector&           random_lower         ,
         const d_vector&           random_upper         ,
         const d_vector&           random_in            ,
         cppad_mixed&              mixed_object
      );
      //
      bool get_nlp_info(
         Index&          n            ,
         Index&          m            ,
         Index&          nnz_jac_g    ,
         Index&          nnz_h_lag    ,
         IndexStyleEnum& index_style
      ) override;
      bool get_bounds_info(
            Index       n        ,
            Number*     x_l      ,
            Number*     x_u      ,
            Index       m        ,
            Number*     g_l      ,
            Number*     g_u
      ) override;
      bool get_starting_point(
         Index           n            ,
         bool            init_x       ,
         Number*         x            ,
         bool            init_z       ,
         Number*         z_L          ,
         Number*         z_U          ,
         Index           m            ,
         bool            init_lambda  ,
         Number*         lambda
      ) override;
      bool eval_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number&         obj_value
      ) override;
      bool eval_grad_f(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Number*         grad_f
      ) override;
      bool eval_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Number*         g
      ) override;
      bool eval_jac_g(
         Index           n        ,
         const Number*   x        ,
         bool            new_x    ,
         Index           m        ,
         Index           nele_jac ,
         Index*          iRow     ,
         Index*          jCol     ,
         Number*         values
      ) override;
      bool eval_h(
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
      ) override;
      void finalize_solution(
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
      ) override;
      bool intermediate_callback(
         AlgorithmMode               mode,
         Index                       iter,
         Number                      obj_value,
         Number                      inf_pr,
         Number                      inf_du,
         Number                      mu,
         Number                      d_norm,
         Number                      regularization_size,
         Number                      alpha_du,
         Number                      alpha_pr,
         Index                       ls_trials,
         const IpoptData*            ip_data,
         IpoptCalculatedQuantities*  ip_cq
      ) override;
      // -----------------------------------------------------------------
      bool adapt_derivative_chk( bool trace, double relative_tol);
      void fixed_eq_constrain(const d_vector& fixed_opt);
   };
} } // END_CPPAD_MIXED_NAMESPACE

# endif
