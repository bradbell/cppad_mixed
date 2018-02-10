// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_IPOPT_FIXED_HPP
# define CPPAD_MIXED_IPOPT_FIXED_HPP

/*
$begin ipopt_fixed$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	Ipopt
	nlp
	inf
	std
$$

$section Ipopt NLP Class Used to Optimize Fixed Effects$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

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

$childtable%src/ipopt_fixed.cpp
	%src/ipopt_fixed/ctor.cpp
	%src/ipopt_fixed/get_nlp_info.cpp
	%src/ipopt_fixed/get_bounds_info.cpp
	%src/ipopt_fixed/get_starting_point.cpp
	%src/ipopt_fixed/eval_f.cpp
	%src/ipopt_fixed/eval_grad_f.cpp
	%src/ipopt_fixed/eval_g.cpp
	%src/ipopt_fixed/eval_jac_g.cpp
	%src/ipopt_fixed/eval_h.cpp
	%src/ipopt_fixed/finalize_solution.cpp
	%example/ipopt_xam.omh
%$$

$end
-----------------------------------------------------------------------------
*/
# include <coin/IpTNLP.hpp>
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
		// ---------------------------------------------------------------
		// member variables set during constructor
		//
		const std::string& random_ipopt_options_;
		const double       fixed_tolerance_;  // ipopt relative tolerance
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
		// random constraint jacobian
		CppAD::mixed::sparse_rcv ran_con_jac_rcv_;
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
		// adaptive_derivative_check
		bool adaptive_called_;
		// The rest of these variables only set by adaptive_derivative_check
		//
		// scale factor (multiplies) components of f(x)
		double scale_f_;
		//
		// scale factor for components of g(x)
		d_vector scale_g_;
		//
		// maps jac_g sparsity index to row and column index in g'(x).
		s_vector jac_g_row_, jac_g_col_;
		// ---------------------------------------------------------------
		// temporaries (size set by constructor only)
		d_vector        fixed_tmp_;         // size n_fixed_
		d_vector        c_vec_tmp_;         // size n_fix_con_
		d_vector        A_uhat_tmp_;        // size n_ran_con_
		d_vector        H_beta_tmp_;        // size n_fixed_
		d_vector        w_fix_con_tmp_;     // size n_fix_con_
		d_vector        w_laplace_obj_tmp_;  // size n_ran_con_ + 1

		// this vector has size fix_likelihood_nabs_ + 1
		d_vector        w_fix_likelihood_tmp_;

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
		//  get and clear the current ipopt_fixed error message
		std::string get_error_message(void) const
		{	return error_message_; }
		//
		// get the fixed effects corresponding to theh error message
		double get_error_fixed(size_t j) const
		{	return error_fixed_[j]; }
		//
		void clear_error_message(void)
		{	error_message_ = ""; }
		//
		// get minus infinity
		double nlp_lower_bound_inf(void) const
		{	return nlp_lower_bound_inf_; }
		//
		// get plus infinity
		double nlp_upper_bound_inf(void) const
		{	return nlp_upper_bound_inf_; }
		//
		// optimal solution
		CppAD::mixed::fixed_solution solution(void) const
		{	return solution_; }
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
		virtual bool get_scaling_parameters(
			Number&            obj_scaling    ,
			bool&              use_x_scaling  ,
			Index              n              ,
			Number*            x_scaling      ,
			bool&              use_g_scaling  ,
			Index              m              ,
			Number*            g_scaling
		);
		// -----------------------------------------------------------------
		bool adaptive_derivative_check( bool trace, double relative_tol);
	};
} } // END_CPPAD_MIXED_NAMESPACE

# endif
