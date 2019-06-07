// ------------------------------------------------------------ 
// Copyright (C) Bradley M. Bell 1998-2015, All rights reserved 
// ------------------------------------------------------------ 
Keyword = 
[
'cppad_mixed  C++ Laplace Approximation of Mixed Effects Models: cppad_mixed-20190607  ',' license source code repository notation fixed theta random data z prior density p(theta) p(z|theta) p(u|theta) p(y|thetau) constraint function c(theta) optimal u^(theta) matrix a*u^(theta) problem maximum likelihood constraints negative log-density vector ',
'install_unix  Installing cppad_mixed in Unix  ',' system requirements c++ compiler git cmake pkg-config wget fortran gsl suitesparse download paths pkg_config_path ld_library_path special run_cmake.sh eigen ipopt command check speed example installation linking using ',
'example_install.sh  An Example Installation  ',' syntax existing source ',
'run_cmake.sh  bin/run_cmake.sh: User Configuration Options  ',' verbose_makefile build_type cppad_prefix eigen_prefix ipopt_prefix debug release extra_cxx_flags cmake_libdir ldlt_cholmod optimize_cppad_function for_hes_sparsity testing speed memory ',
'check_install.sh  Example and Test Using the Installed Version of cppad_mixed  ',' syntax build_type prefixes cmake_libdir example_file pkg_config_path ld_library_path create temporary example_name main gsl_libs ipopt_libs suitesparse_libs compile link run ',
'theory  Laplace Approximation for Mixed Effects Models  ',' reference total likelihood random f(theta u) assumption fixed g(theta) optimal u^(theta) objective h(theta r(theta) l(theta) derivative constraints approximate first order u(beta second w(beta h(beta function b(beta hessian sparse observed information ',
'base_class  The cppad_mixed Base Class  ',' public private ',
'public  cppad_mixed: Public Declarations  ',' types user defined functions ran_likelihood fix_likelihood fix_constraint fatal_error warning constructor destructor initialize optimize_random optimize_fixed information_mat sample_fixed sample_random ',
'private  cppad_mixed: Private Declarations  ',' n_fixed_ n_random_ quasi_fixed_ bool_sparsity_ a_rcv_ initialize_done_ cppad_error_handler ran_like_fun_ ran_jac_fun_ ran_hes_fun_ ldlt_ran_hes_ hes_cross_ laplace_obj_fun_ laplace_obj_hes_ fix_like_fun_ fix_like_jac_ fix_like_hes_ fix_con_fun_ fix_con_jac_ fix_con_hes_ template member functions pack unpack initialization init_ldlt_ran_hes init_fix_con init_fix_like init_hes_cross init_ran_jac init_ran_hes init_laplace_obj_hes init_ran_like try try_initialize try_optimize_random try_optimize_fixed try_information_mat try_sample_fixed sample_random other fix_con_eval fix_like_eval logdet_jac ran_con_eval ran_con_jac ran_obj_eval ran_obj_jac update_factor ',
'derived_ctor  User Defined Class Derived From cppad_mixed  ',' syntax see also mixed_derived mixed_object n_fixed n_random quasi_fixed true false bool_sparsity a_rcv ... errorhandler example ',
'derived_ctor.cpp  mixed_cppad Derived Class: Example and Test  ',' ',
'ran_likelihood  User Defined Random Likelihood Function  ',' syntax public mixed_object virtual fixed_vec random_vec constant default example ',
'ran_likelihood.cpp  Random Likelihood: Example and Test  ',' ',
'fix_likelihood  User Defined Fixed Likelihood Function  ',' syntax public mixed_object a1_double virtual fixed_vec constant default example ',
'fix_likelihood.cpp  Random Likelihood: Example and Test  ',' ',
'fix_constraint  User Defined Fixed Effects Constraint Function  ',' syntax public mixed_object a1_double virtual fixed_vec default example ',
'fix_constraint.cpp  Using Constraints: Example and Test  ',' model ',
'initialize  Initialization After Constructor  ',' syntax public purpose mixed_object fixed_vec random_vec size_map init_newton_checkpoint_done_ example ',
'optimize_random  Optimize Random Effects  ',' syntax public purpose mixed_object options evaluation_method fixed_vec random_lower random_upper random_in random_out example ',
'optimize_random.cpp  Optimize Random Effects: Example and Test  ',' ',
'optimize_fixed  Optimize Fixed Effects  ',' syntax public purpose inf mixed_object fixed_ipopt_options derivative_test hessian_approximation max_iter accept_after_max_steps nlp_scaling_method random_ipopt_options fixed_lower fixed_upper fix_constraint_lower fix_constraint_upper fixed_scale fixed_in random_lower random_upper random_in solution laplace example ipopt_fixed ',
'optimize_fixed.cpp  Optimize Fixed Effects: Example and Test  ',' model objective first order partials source code ',
'ipopt_options  An Ipopt Options Argument  ',' prototype format string integer numeric ',
'ipopt_trace  Description of Ipopt Tracing Output  ',' iter objective inf_pr inf_du lg(mu) ||d|| lg(rg) alpha_du alpha_pr ls reference ',
'information_mat  Compute the Observed Information For Fixed Effects  ',' syntax purpose quasi_fixed mixed_object solution random_opt information_rcv example ',
'information_mat.cpp  Observed Information Matrix: Example and Test  ',' model ',
'sample_fixed  Sample Posterior for Fixed Effects  ',' syntax see also prototype public purpose manage_gsl_rng mixed_object information_rcv solution fixed_lower fixed_upper random_opt positive definite theory notation unconstrained subset covariance approximate constraint equations implicit example other method 2do ',
'sample_fixed.cpp  Sample From Fixed Effects Posterior: Example and Test  ',' ',
'sample_conditional  Sample Posterior for Fixed Effects Using Conditional Covariance  ',' syntax replaced prototype public purpose manage_gsl_rng mixed_object information_info solution fixed_lower fixed_upper random_opt theory notation subset unconstrained constraint equations example ',
'sample_random  Simulation the Posterior Distribution for Random Effects  ',' syntax see also prototype public purpose manage_gsl_rng mixed_object random_ipopt_options fixed_vec random_lower random_upper random_in covariance example ',
'sample_random.cpp  Sample From Fixed Effects Posterior: Example and Test  ',' ',
'pack  Pack Fixed Effect and Random Effects Into One Vector  ',' syntax private mixed_object float_unpack float_pack fixed_one fixed_two random_vec both_vec three_vec ',
'unpack  Pack Fixed Effect and Random Effects Into One Vector  ',' syntax private mixed_object float_pack float_unpack fixed_one fixed_two random_vec both_vec three_vec ',
'init_ran_jac  Initialize Jacobian of Random Likelihood w.r.t. Random Effects  ',' syntax private assumptions init_ran_jac_done_ mixed_object fixed_vec random_vec ran_jac_a1fun_ ran_jac2hes_rc_ ',
'ran_jac_fun.cpp  ran_jac_fun: Example and Test  ',' private ',
'init_ran_hes  Initialize Hessian of Random Likelihood w.r.t Random Effects  ',' syntax private assumptions init_ran_hes_done_ mixed_object fixed_vec random_vec ran_hes_uu_rcv_ index order ran_hes_fun_ ',
'ran_hes_fun.cpp  ran_hes_fun_: Example and Test  ',' private ',
'init_laplace_obj  Second Order Representation of Laplace Objective and Constraints  ',' syntax private assumptions init_laplace_obj_done_ mixed_object fixed_vec random_vec laplace_obj_fun_ beta theta range space ',
'init_ldlt_ran_hes  Initialize Cholesky Factor of Hessian of Random Likelihood  ',' syntax private assumptions init_ldlt_ran_hes_done_ mixed_object a1_ldlt_ran_hes_ ',
'init_fix_con  Initialize Constraints as Function of Fixed Effects  ',' syntax private init_fix_con_done_ mixed_object fixed_vec fix_con_fun_ fix_con_jac_ fix_con_hes_ ',
'init_fix_like  Initialize Fixed Likelihood  ',' syntax private init_fix_like_done_ mixed_object fixed_vec fix_like_fun_ fix_like_jac_ fix_like_hes_ ',
'init_hes_cross  Cross Terms of Sparse Hessian w.r.t Fixed and Random Effects  ',' syntax private assumptions init_hes_cross_done_ mixed_object fixed_vec random_vec ran_like_fun_ ran_like_a1fun_ order example ',
'hes_cross.cpp  Hessian Cross Terms: Example and Test  ',' private ',
'init_laplace_obj_hes  Initialize Hessian of Approximate Laplace Objective  ',' syntax private assumptions init_laplace_obj_hes_done_ mixed_object fixed_vec random_vec laplace_obj_fun_ ',
'init_ran_like  Initialize Random Likelihood  ',' syntax private init_ran_like_done_ mixed_object fixed_vec random_vec ran_like_fun_ ran_like_a1fun_ ',
'fix_con_eval  Evaluate Fixed Constraint Function  ',' syntax private mixed_object fixed_vec example ',
'fix_con_eval.cpp  constraint_eval: Example and Test  ',' private ',
'fix_con_hes  Hessian of Fixed Constraints  ',' syntax private mixed_object fixed_vec weight row_out col_out val_out example ',
'fix_con_hes.cpp  constraint_hes: Example and Test  ',' private ',
'fix_con_jac  Jacobian of Fixed Constraint  ',' syntax private mixed_object fixed_vec row_out col_out val_out example ',
'fix_con_jac.cpp  constraint_jac: Example and Test  ',' private ',
'fix_like_eval  Evaluate Fixed Likelihood  ',' syntax private mixed_object fix_likelihood fixed_vec example ',
'fix_like_eval.cpp  fix_like_eval: Example and Test  ',' private ',
'fix_like_hes  Hessian of Fixed Likelihood  ',' syntax private mixed_object fix_likelihood fixed_vec weight row_out col_out val_out example ',
'fix_like_hes.cpp  fix_like_hes: Example and Test  ',' private ',
'fix_like_jac  Jacobian of Fixed Likelihood  ',' syntax private mixed_object fix_likelihood fixed_vec row_out col_out val_out example ',
'fix_like_jac.cpp  fix_like_jac: Example and Test  ',' private ',
'logdet_jac  Jacobian of Log Determinant of Hessian w.r.t. Random Effects  ',' syntax private purpose mixed_object ldlt_ran_hes_ fixed_vec random_vec logdet_fix logdet_ran example ',
'logdet_jac.cpp  logdet_jac: Example and Test  ',' private ',
'ran_like_hes  Hessian of Random Likelihood w.r.t. Random Effects  ',' syntax prototype private assumptions purpose n_fixed n_random ran_jac_a1fun ran_hes_uu_rc theta_u ran_hes_uu_rcv example ',
'ran_like_hes.cpp  ran_like_hes: Example and Test  ',' private ',
'ran_con_eval  Evaluate the Random Constraint Function  ',' syntax private mixed_object random_vec au example ',
'ran_con_eval.cpp  ran_con_eval: Example and Test  ',' private model ',
'ran_con_jac  Jacobian of the Random Constraint Function  ',' syntax private mixed_object fixed_vec random_vec jac_rcv sparsity pattern sparse matrix column major order example ',
'ran_con_jac.cpp  ran_con_jac: Example and Test  ',' private model ',
'ran_obj_eval  Evaluate Laplace Approximation and Laplace Objective  ',' syntax private purpose mixed_object ldlt_ran_hes_ fixed_vec random_vec example ',
'ran_obj_eval.cpp  ran_obj_eval: Example and Test  ',' private model ',
'ran_obj_jac  Derivative of Laplace Objective  ',' syntax private purpose mixed_object ldlt_ran_hes_ fixed_vec random_vec r_fixed example ',
'ran_obj_jac.cpp  ran_obj_jac: Example and Test  ',' private model ',
'laplace_obj_hes  Hessian of Laplace Objective and Random Constraints  ',' syntax private purpose mixed_object fixed_vec random_vec weight row_out col_out val_out example ',
'laplace_obj_hes.cpp  laplace_obj_hes: Example and Test  ',' private model ',
'update_factor  Update the Factorization of Hessian w.r.t. Random Effects  ',' syntax private purpose mixed_object fixed_vec random_vec ran_hes_fun_ ldlt_ran_hes_ example ',
'update_factor.cpp  update_factor: Example and Test  ',' private ',
'namespace  The CppAD::mixed Namespace Declarations  ',' configure.hpp public private ',
'typedef  Types Defined in the CppAD Mixed Namespace  ',' syntax begin scalar a1_double vector s_vector d_vector a1_vector sparse sparse_rc d_sparse_rcv a1_sparse_rcv end ',
'configure.hpp  Preprocessor Configuration Symbols in configure.hpp  ',' private cppad_mixed_version cppad_mixed_for_hes_sparsity cppad_mixed_optimize_ad_function cppad_mixed_ldlt_cholmod cppad_mixed_ldlt_class cppad_mixed_null_ptr ',
'exception  CppAD Mixed Exceptions  ',' syntax public thrower brief catcher description ',
'fixed_solution  Optimal Solution Returned by optimize_fixed  ',' syntax public convention fixed_opt fixed_lag bounds fix_con_lag ran_con_lag ',
'ipopt_app_status  Map Ipopt Application Return Status to a String  ',' syntax public prototype ',
'ipopt_fixed  Ipopt NLP Class Used to Optimize Fixed Effects  ',' private get_error_message() clear_error_message() nlp_lower_bound_inf() nlp_upper_bound_inf() solution() default destructor ',
'ipopt_fixed_ctor  Ipopt Fixed Optimization Callback Constructor and Destructor  ',' syntax private references abort_on_eval_error random_ipopt_options fixed_tolerance fixed_lower fixed_upper fix_constraint_lower fix_constraint_upper fixed_scale fixed_in random_lower random_upper random_in mixed_object argument constants n_fixed_ n_random_ n_fix_con_ n_ran_con temporaries fixed_tmp_ c_vec_tmp_ a_uhat_tmp_ h_beta_tmp_ w_fix_con_tmp_ w_laplace_obj_tmp_ w_fix_likelihood_tmp_ fix_likelihood_vec_tmp_ effectively member variables nlp_lower_bound_inf_ nlp_upper_bound_inf_ fix_likelihood_nabs_ nnz_jac_g_ nnz_h_lag_ lag_hes_row_ lag_hes_col_ laplace_obj_hes_2_lag_ fix_like_hes_2_lag_ sparsity information ran_con_jac_rcv_ laplace_obj_hes_info_ adaptive_called_ prototype ',
'ipopt_fixed_get_nlp_info  Return Information About Problem Sizes  ',' syntax nnz_jac_g nnz_h_lag index_style ok prototype ',
'ipopt_fixed_get_bounds_info  Return Optimization Bounds  ',' syntax x_l x_u g_l g_u ok prototype ',
'ipopt_fixed_get_starting_point  Return Initial Values Where Optimization is Started  ',' syntax init_x init_z z_l z_u init_lambda ok prototype ',
'ipopt_fixed_eval_f  Compute Value of Objective  ',' syntax new_x obj_val ok prototype ',
'ipopt_fixed_eval_grad_f  Compute Gradient of the Objective  ',' syntax new_x ok prototype ',
'ipopt_fixed_eval_g  Compute Value of Constraint Functions  ',' syntax new_x ok prototype ',
'ipopt_fixed_eval_jac_g  Compute Jacobian of Constraint Functions  ',' syntax new_x nele_jac irow jcol values ok prototype ',
'ipopt_fixed_eval_h  Compute the Hessian of the Lagrangian  ',' syntax mixed_object.quasi_fixed_ new_x obj_factor lambda new_lambda nele_hess irow jcol values ok prototype ',
'ipopt_fixed_finalize_solution  Get Solution Results  ',' syntax solution_ z_l z_u m lambda obj_value ip_data ip_cq status success maxiter_exceeded cputime_exceeded stop_at_tiny_step stop_at_acceptable_point local_infeasibility user_requested_stop diverging_iterates restoration_failure error_in_step_computation invalid_number_detected prototype ',
'ipopt_fixed_adaptive_derivative_check  Adaptive Step Size check of Derivative Calculations  ',' syntax trace relative_step relative_tol infinity ok error_message_ adaptive_called_ fixed_scale_ scale_f_ scale_g_ jac_g_row_ jac_g_col_ prototype ',
'ipopt_fixed_new_random  Compute New Random Effects and Update Factor  ',' syntax n_random_ random_ipopt_options_ random_lower_ random_upper_ random_in_ fixed_vec random_cur_ mixed_object_ prototype ',
'ipopt_xam  Example Use of Ipopt  ',' problem lagrangian stationary conditions solution ',
'ipopt_nlp_xam  Ipopt Example: Declare Non-linear Program Problem Class  ',' ',
'ipopt_xam_ctor  Ipopt Example: Constructor and Destructor  ',' ',
'ipopt_xam_get_nlp_info  Return Information About Problem Sizes  ',' syntax nnz_jac_g nnz_h_lag index_style ok source ',
'ipopt_xam_get_bounds_info  Return Optimization Bounds  ',' syntax x_l x_u g_l g_u ok source ',
'ipopt_xam_get_starting_point  Return Initial Values Where Optimization is Started  ',' syntax init_x init_z z_l z_u init_lambda ok source ',
'ipopt_xam_eval_f  Compute Value of Objective  ',' syntax new_x obj_val ok source ',
'ipopt_xam_eval_grad_f  Compute Gradient of the Objective  ',' syntax new_x ok source ',
'ipopt_xam_eval_g  Compute Value of Constraint Functions  ',' syntax new_x ok source ',
'ipopt_xam_eval_jac_g  Compute Jacobian of Constraint Functions  ',' syntax new_x nele_jac irow jcol values ok source ',
'ipopt_xam_eval_h  Compute the Hessian of the Lagrangian  ',' syntax new_x obj_factor lambda new_lambda nele_hess irow jcol values ok source ',
'ipopt_xam_finalize_solution  Get Solution Results  ',' syntax z_l z_u lambda obj_value ip_data ip_cq status success maxiter_exceeded cputime_exceeded stop_at_tiny_step stop_at_acceptable_point local_infeasibility user_requested_stop diverging_iterates restoration_failure error_in_step_computation invalid_number_detected source ',
'ipopt_xam_intermediate_callback  Ipopt Example: Optimization Progress Report  ',' syntax mode iter obj_value inf_pr inf_du mu d_norm regularization_size alpha_du alpha_pr ls_trials ok source ',
'ipopt_run_xam  Ipopt: Example and Test  ',' syntax ok source ',
'ipopt_random  Ipopt NLP Class Used to Optimize Random Effects  ',' syntax private get_error_message() clear_error_message() nlp_lower_bound_inf() nlp_upper_bound_inf() solution() default destructor ',
'ipopt_random_ctor  Ipopt Random Optimization Callback Constructor  ',' syntax fixed_vec random_lower random_upper random_in mixed_object member variables nlp_lower_bound_inf_ nlp_upper_bound_inf_ nnz_h_lag_ ',
'ipopt_random_get_nlp_info  Return Information About Problem Sizes  ',' syntax nnz_jac_g nnz_h_lag index_style ok prototype ',
'ipopt_random_get_bounds_info  Return Optimization Bounds  ',' syntax x_l x_u g_l g_u ok prototype ',
'ipopt_random_get_starting_point  Return Initial Values Where Optimization is Started  ',' syntax init_x init_z z_l z_u init_lambda ok prototype ',
'ipopt_random_eval_f  Compute Value of Objective  ',' syntax new_x obj_val ok objective_current_ mixed_object_.ran_like_fun_ prototype ',
'ipopt_random_eval_grad_f  Compute Gradient of the Objective  ',' syntax new_x ok mixed_object_.ran_like_fun_ prototype ',
'ipopt_random_eval_g  Compute Value of Constraint Functions  ',' syntax new_x ok mixed_object_.ran_like_fun_ prototype ',
'ipopt_random_eval_jac_g  Compute Jacobian of Constraint Functions  ',' syntax new_x nele_jac irow jcol values ok mixed_object_.ran_like_fun_ prototype ',
'ipopt_random_eval_h  Compute the Hessian of the Lagrangian  ',' syntax mixed_object.quasi_fixed_ new_x obj_factor lambda new_lambda nele_hess irow jcol values ok mixed_object_.ran_like_fun_ mixed_object_.ran_hes_fun_ prototype ',
'ipopt_random_finalize_solution  Get Solution Results  ',' syntax solution_ z_l z_u lambda obj_value ip_data ip_cq status success maxiter_exceeded cputime_exceeded stop_at_tiny_step stop_at_acceptable_point local_infeasibility user_requested_stop diverging_iterates restoration_failure error_in_step_computation invalid_number_detected prototype ',
'ldlt_cholmod  A Cholmod Cholesky Factor Class  ',' see also private purpose factorization example preprocessor symbols ',
'ldlt_cholmod_ctor  Cholmod LDLT Constructor  ',' syntax prototype private nrow_ pointers common_ simplicial factorization ldl\' example ',
'ldlt_cholmod_dtor  Cholmod Destructor  ',' private discussion ',
'ldlt_cholmod_init  Initialize Factor for a Specific Sparsity Pattern  ',' syntax prototype private ldlt_obj h_rc h_rc_ assumptions sym_matrix_ factor_ rhs_ rhs_set_ order operations example ',
'ldlt_cholmod_pattern  Update Factorization Using new Matrix Values  ',' syntax prototype private ldlt_obj h_rc order operations example ',
'ldlt_cholmod_update  Update Factorization Using new Matrix Values  ',' syntax prototype private purpose ldlt_obj h_rcv sym_matrix_ factor_ ok order operations example ',
'ldlt_cholmod_logdet  Compute Log Determinant for Current Factor  ',' syntax prototype private ldlt_obj negative example ',
'ldlt_cholmod_solve_H  Solve Linear Equations Using Factor  ',' syntax prototype private purpose ldlt_obj row val_in val_out example ',
'ldlt_cholmod_sim_cov  Simulations with Covariance Corresponding to Factored Matrix  ',' syntax prototype private purpose ldlt_obj positive definite ok example ',
'ldlt_cholmod_inv  Compute a Subset of the Inverse of Factored Matrix  ',' under construction syntax prototype private purpose ldlt_obj row_in col_in val_out method member variables sparseinv_p_ sparseinv_i_ out2sparseinv_order_ ',
'ldlt_cholmod.cpp  Example Using Cholmod LDLT Factorization  ',' problem description constructor init pattern update logdet solve_h sim_cov source code ',
'cholmod_solve_xam  Example Using cholmod_solve With a Non-positive Matrix  ',' problem description source code ',
'cholmod_solve2_a.cpp  Example Using cholmod_solve2 Cholesky Factorization  ',' problem description bset source code ',
'cholmod_solve2_sim.cpp  Cholmod Posterior Simulations Using Sparse Hessian of Likelihood  ',' problem example source code ',
'ldlt_eigen  An Eigen LDLT Factor Class  ',' see also private factorization h double example ',
'ldlt_eigen_ctor  Eigen LDLT Constructor  ',' syntax prototype private n_row_ ptr_ ',
'ldlt_eigen_init  Initialize LDLT Factor for a Specific Sparsity Pattern  ',' syntax prototype private ldlt_obj h_rc h_rc_ order operations example ',
'ldlt_eigen_pattern  Update Factorization Using new Matrix Values  ',' syntax prototype private ldlt_obj h_rc order operations example ',
'ldlt_eigen_update  Update Factorization Using new Matrix Values  ',' syntax prototype private purpose ldlt_obj h_rcv ptr_ ok order operations example ',
'ldlt_eigen_split  Split Out a Factorization  ',' syntax prototype private purpose ldlt_obj example ',
'ldlt_eigen_logdet  Compute Log Determinant for Current LDLT Factor  ',' syntax prototype private ldlt_obj negative example ',
'ldlt_eigen_solve_H  Solve Linear Equations Using Stored Factor  ',' syntax prototype private purpose ldlt_obj row val_in val_out example ',
'ldlt_eigen_sim_cov  Simulations with Covariance Corresponding to Factored Matrix  ',' syntax prototype private purpose ldlt_obj positive definite ok example ',
'ldlt_eigen_inv  Compute a Subset of the Inverse of Factored Matrix  ',' syntax prototype private purpose ldlt_obj row_in col_in val_out method ',
'ldlt_eigen_solve_LDLT  Solve Linear Equations Corresponding to L, D, and P Factors  ',' syntax prototype private purpose b example ',
'ldlt_eigen.cpp  Example Using Eigen LDLT Factorization  ',' problem description constructor init pattern update logdet solve_h split solve_ldlt sim_cov source code ',
'manage_gsl_rng  Set, Get, And Free A GSL Random Number Generator  ',' syntax public purpose new_gsl_rng s_in s_out get_gsl_rng free_gsl_rng example ',
'manage_gsl_rng.cpp  Manage GSL Random Number Generator: Example and Test  ',' ',
'sparse_hes_rcv  Sparse Hessian Computation Structure  ',' syntax private purpose subset nnz row col val work computing hessians f not_used_pattern not_used_coloring ',
'sparse_hes_info  Sparse Hessian Information  ',' syntax private purpose row col val work call not_used hes_info.val ',
'sparse_jac_rcv  Sparse Jacobian Computation Structure  ',' syntax private purpose subset nnz row col val forward work computing jacobians mode reverse group_max not_used_pattern not_used_coloring ',
'sparse_mat_info  Sparse Matrix Information  ',' syntax public purpose row k col val resize notation sparsity pattern empty column major order lower triangular ',
'triple2eigen  Convert Row, Column, Value Triple to an Eigen Sparse Matrix  ',' syntax prototype private scalar nr nc sparsity pattern ',
'order2random  Second Order Representation of Random Effects, W(beta, theta, u)  ',' syntax prototype private n_fixed n_random jac_a1fun a1_ldlt_ran_hes theta_u example ',
'order2random.cpp  order2random: Example and Test  ',' private ',
'undetermined  Express An Undetermined Linear System As Dependent and Independent Variables  ',' syntax prototype private purpose nr nc tol rank example 2do ',
'undetermined.cpp  undetermined: Example and Test  ',' private ',
'sparse_low_tri_sol  Solve a Sparse Lower Triangular Linear System  ',' syntax prototype private left right result example ',
'sparse_low_tri_sol.cpp  sparse_low_tri_sol: Example and Test  ',' private problem ',
'sparse_up_tri_sol  Solve a Sparse Upper Triangular Linear System  ',' syntax prototype private left right result example ',
'sparse_up_tri_sol.cpp  sparse_up_tri_sol: Example and Test  ',' private problem ',
'sparse_scale_diag  Scales the Diagonal of an Eigen Sparse Matrix  ',' syntax purpose private scalar options index example ',
'sparse_scale_diag.cpp  sparse_scale_diag: Example and Test  ',' private description ',
'sparse_low2sym  Convert an Eigen Lower Triangular Matrix To a Symmetric Matrix  ',' syntax private scalar options index example ',
'sparse_low2sym.cpp  sparse_low2sym: Example and Test  ',' private ',
'sparse_mat2low  Extract the Lower Triangular From an Eigen Symmetric Matrix  ',' syntax private scalar options index example ',
'sparse_mat2low.cpp  sparse_mat2low: Example and Test  ',' private ',
'sparse_eigen2info  Convert An Eigen Sparse Matrix to a sparse_mat_info Representation  ',' syntax private option index empty input non-empty example ',
'sparse_eigen2info.cpp  sparse_eigen2info: Example and Test  ',' private ',
'sparse_info2eigen  Convert a sparse_mat_info Representation to An Eigen Sparse Matrix  ',' syntax private option index example ',
'sparse_info2eigen.cpp  sparse_info2eigen: Example and Test  ',' private ',
'sparse_eigen2rcv  Convert an Eigen Sparse Matrix to a Cppad Sparse Matrix  ',' syntax prototype private scalar m_eigen m_rcv example ',
'sparse_eigen2rcv.cpp  sparse_eigen2rcv: Example and Test  ',' private ',
'sparse_rcv2eigen  Convert a CppAD Sparse Matrix to an Eigen Sparse Matrix  ',' syntax prototype private scalar m_rcv m_eigen example ',
'sparse_rcv2eigen.cpp  sparse_rcv2eigen: Example and Test  ',' private ',
'sparse_print  Print and Eigen Sparse Matrix  ',' syntax private scalar label ',
'sparsity_print  Print a CppAD Internal Sparsity Pattern  ',' syntax private label ',
'sparse_ad_cholesky  Sparse Cholesky Factorization as an Atomic CppAD Operation  ',' syntax not used purpose notation private public member functions type declarations variables ',
'sparse_ad_cholesky_initialize  Initialize Sparse AD Cholesky Factorization  ',' syntax public / private ad_alow restriction ',
'sparse_ad_cholesky_p  Using Sparse AD Cholesky Permutation P  ',' syntax prototype public / private example ',
'sparse_ad_cholesky_eval  Using Sparse AD Cholesky Factor L  ',' syntax public / private ad_alow ad_l example ',
'set_jac_sparsity  Set the Jacobian Sparsity Pattern  ',' syntax private n_set end elements ',
'set_hes_sparsity  Set the Hessian Sparsity Pattern  ',' syntax purpose private jac_sparsity n_set end elements ',
'sparse_ad_chol_eval.cpp  Sparse AD Cholesky Factorization: Example and Test  ',' problem source ',
'sparse_ad_chol_perm.cpp  Sparse AD Cholesky Permutation: Example and Test  ',' source ',
'sparse_ad_chol_eq.cpp  Using Sparse AD Cholesky To Solve Equations: Example and Test  ',' source ',
'sparse_ad_chol_var.cpp  Sparse AD Cholesky Variable Calculation: Example and Test  ',' problem permutation factor source ',
'sparse_ad_chol_sp1.cpp  Sparse AD Cholesky Sparsity Calculation: Example and Test  ',' problem permutation factor source ',
'sparse_ad_chol_sp2.cpp  Sparse AD Cholesky Sparsity Calculation: Example and Test  ',' problem permutation factor source ',
'user  User API Examples  ',' functions defined called cppad_mixed programs that are speed tests demonstrate specific features ',
'speed  Example Programs That are Also Speed and Memory Tests  ',' purpose ',
'ar1_xam.cpp  A First Order Auto-Regressive Example and Speed Test  ',' syntax problem data p( y_t | theta ) command arguments random_seed number_random quasi_fixed trace_optimize_fixed ipopt_solve bool_sparsity hold_memory derivative_test start_near_solution output cppad_mixed_version ldlt_cholmod optimize_cppad_function actual_seed initialize_bytes initialize_seconds optimize_fixed_seconds optimize_random_seconds information_mat_seconds sample_fixed_seconds final_bytes theta_0_estimate ar1_xam_ok source code ',
'ar1_xam.sh  Example Using ar1_xam  ',' syntax test2run normal callgrind massif source code ',
'capture_xam.cpp  A Capture Example and Speed Test  ',' syntax reference command arguments random_seed number_random quasi_fixed trace_optimize_fixed ipopt_solve bool_sparsity hold_memory derivative_test start_near_solution number_fixed_samples number_locations max_population mean_population mean_logit_probability std_logit_probability random_constraint output cppad_mixed_version ldlt_cholmod optimize_cppad_function actual_seed initialize_bytes initialize_seconds optimize_fixed_seconds optimize_random_seconds information_mat_seconds sample_fixed_seconds final_bytes sum_random_effects mean_population_estimate mean_logit_probability_estimate std_logit_probability_estimate mean_population_std mean_logit_probability_std std_logit_probability_std mean_population_ratio mean_logit_probability_ratio std_probability_ratio capture_xam_ok notation p(y_it|n_iq_t) p(n_i|theta) q_t(thetau) m_i p(y_i|thetau) p(y|thetau) p(u|theta) p(theta) p(z|theta) c(theta) source code ',
'capture_xam.sh  Example Using capture_xam  ',' syntax test2run normal callgrind massif source code ',
'abs_density.cpp  Absolute Value In Log-Density: Example and Test  ',' model ',
'no_random.cpp  No Random Effects: Example and Test  ',' model ',
'ran_constraint.cpp  Constraints On Random Effects: Example and Test  ',' ',
'lasso.cpp  Lasso on Fixed Effects: Example and Test  ',' model ',
'data_mismatch.cpp  Random Effects Variance May Cause Data Mismatch  ',' model theory derivatives objective ',
'opt_ran_nan.cpp  Nan\'s During Optimization of Random Effects: Example and Test  ',' ',
'whats_new  Changes and Additions to cppad_mixed  ',' ',
'whats_new_19  Changes and Additions to cppad_mixed During 2019  ',' 06 ',
'whats_new_18  Changes and Additions to cppad_mixed During 2018  ',' 10-08 09-25 08-30 08-27 08-22 08-20 08-18 08-09 08-08 07-25 07-12 06-30 06-29 06-20 06-14 06-09 06-04 05-21 05-07 05-03 04-06 03-22 03-10 02-20 02-12 02-11 02-10 02-08 02-07 02-05 01-23 01-22 01-21 01-15 01-14 ',
'whats_new_17  Changes and Additions to cppad_mixed During 2017  ',' 12-28 12-25 12-22 12-16 12-10 10-27 10-24 10-09 10-07 09-30 09-23 09-21 09-18 09-16 09-15 09-14 09-02 08-30 08-01 04-24 04-23 04-06 04-02 03-27 03-25 03-23 03-20 03-12 03-11 03-10 03-09 03-08 03-06 03-02 03-01 01-26 01-24 01-22 01-14 ',
'whats_new_16  Changes and Additions to cppad_mixed During 2016  ',' 11-09 11-07 11-04 11-02 10-30 10-28 10-27 10-25 10-24 10-18 10-16 10-15 10-14 10-06 09-30 09-27 09-24 09-23 07-28 07-27 07-26 07-25 07-20 07-18 07-14 07-13 07-12 07-10 07-09 06-24 06-22 06-19 06-18 06-17 06-13 06-12 06-11 06-07 06-06 06-05 06-04 06-03 optimize_random optimize_fixed sample_random 05-15 05-11 05-08 05-06 05-04 05-03 04-29 04-27 04-23 04-19 04-18 04-17 04-16 04-10 04-15 04-09 04-08 04-07 04-06 04-05 04-03 04-02 04-01 03-29 03-28 03-09 02-26 02-06 01-26 01-25 01-22 01-21 01-19 01-16 01-15 01-14 01-13 01-10 01-09 01-05 01-04 01-01 ',
'whats_new_15  Changes and Additions to cppad_mixed During 2015  ',' 12-25 12-24 12-16 12-14 12-13 12-10 ',
'wish_list  CppAD Mixed Wish List  ',' multi-threading sparse matrix d_sparse_rcv aborting optimization random constraints windows install ',
'math_notation  Mathematical Notation  ',' b c_l c_u f g p r u^(theta) w y z '
]

var MaxList = 100;
var Nstring = -1;
var Nkeyword = Keyword.length / 2;
Initialize();

function Initialize()
{
	UpdateList();
	document.search.keywords.focus();
}
function UpdateList(event)
{
	key = 0;
	if( window.event )
		key = window.event.keyCode;
	else if( event )
		key = event.which;
	if( key == 13 )
	{	Goto();
		return;
	}
	var string  = document.search.keywords.value;
	if( Nstring == string.length )
		return;
	Nstring     = string.length;

	var word    = string.match(/\S+/g);
	var nword   = 0;
	if(word != null )
		nword   = word.length;

	var pattern = new Array(nword);
	for(var j = 0; j < nword; j++)
		pattern[j] = new RegExp(word[j], 'i');

	var nlist = 0;
	var list  = '';
	for(i = 0; (i < Nkeyword) && (nlist < MaxList) ; i++)
	{
		var match = true;
		for(j = 0; j < nword; j++)
		{	var flag = pattern[j].test(Keyword[2*i]);
			flag     = flag || pattern[j].test(Keyword[2*i+1]);
			match    = match && flag;
		}

		if( match )
		{
			line  = Keyword[2*i].split(/\s+/);
			line  = line.join(' ');
			list  = list + line + '\n';
			nlist = nlist + 1;
		}
	}
	document.search.list.value    = list;
}
function Choose(textarea)
{	var start_select = textarea.value.substring(0, textarea.selectionStart);
	var start_pos    = Math.max(0, start_select.lastIndexOf('\n') );
	var length       = textarea.value.length;
	var end_select   = 
		textarea.value.substring(textarea.selectionEnd, length);
	var end_pos      = end_select.indexOf('\n');
	if( end_pos >= 0 ) 
	{	end_pos = textarea.selectionEnd + end_pos;
	} else 
	{	end_pos = length;
	}
	// highlight the selected line
	textarea.selectionStart = start_pos;
	textarea.selectionEnd   = end_pos;
	// get the choice from the beginning of the line
	var line = textarea.value.substring(start_pos, length);
	var end_choice = line.indexOf(' ');
	if( end_choice >= 0 )
	{	var choice = line.substring(0, end_choice);
		document.search.choice.value = choice.toLowerCase();
	}
	
	return true;
}
function Goto()
{
parent.location = document.search.choice.value + '.htm';
}
