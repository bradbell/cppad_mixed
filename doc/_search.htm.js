// ------------------------------------------------------------ 
// Copyright (C) Bradley M. Bell 1998-2015, All rights reserved 
// ------------------------------------------------------------ 
Keyword = 
[
'cppad_mixed  C++ Laplace Approximation of Mixed Effects Models: cppad_mixed-20190624  ',' license source code repository notation fixed theta random data z prior density p(theta) p(z|theta) p(u|theta) p(y|thetau) constraint function c(theta) optimal u^(theta) matrix a*u^(theta) problem maximum likelihood constraints negative log-density vector ',
'install_unix  Installing cppad_mixed in Unix  ',' system requirements c++ compiler git cmake pkg-config wget fortran gsl suitesparse download paths pkg_config_path ld_library_path special run_cmake.sh eigen ipopt command check speed example installation linking using ',
'example_install.sh  An Example Installation  ',' syntax existing source ',
'run_cmake.sh  bin/run_cmake.sh: User Configuration Options  ',' verbose_makefile build_type cppad_prefix eigen_prefix ipopt_prefix debug release extra_cxx_flags cmake_libdir ldlt_cholmod optimize_cppad_function for_hes_sparsity testing speed memory ',
'check_install.sh  Example and Test Using the Installed Version of cppad_mixed  ',' syntax build_type prefixes cmake_libdir example_file pkg_config_path ld_library_path create temporary example_name main gsl_libs ipopt_libs suitesparse_libs compile link run ',
'theory  Laplace Approximation for Mixed Effects Models  ',' reference total likelihood random f(theta u) assumption fixed g(theta) optimal u^(theta) objective h(theta r(theta) l(theta) derivative constraints approximate first order u(beta second w(beta h(beta function b(beta hessian sparse observed information ',
'base_class  cppad_mixed: Public Declarations  ',' types user defined functions ran_likelihood fix_likelihood fix_constraint fatal_error warning constructor destructor initialize optimize_random optimize_fixed information_mat sample_fixed sample_random ',
'derived_ctor  User Defined Class Derived From cppad_mixed  ',' syntax see also mixed_derived mixed_object n_fixed n_random quasi_fixed true false bool_sparsity a_rcv ... errorhandler example ',
'derived_ctor.cpp  mixed_cppad Derived Class: Example and Test  ',' ',
'ran_likelihood  User Defined Random Likelihood Function  ',' syntax public mixed_object virtual fixed_vec random_vec constant default example ',
'ran_likelihood.cpp  Random Likelihood: Example and Test  ',' ',
'fix_likelihood  User Defined Fixed Likelihood Function  ',' syntax public mixed_object a1_double virtual fixed_vec constant default example ',
'fix_likelihood.cpp  Random Likelihood: Example and Test  ',' ',
'fix_constraint  User Defined Fixed Effects Constraint Function  ',' syntax public mixed_object a1_double virtual fixed_vec default example ',
'fix_constraint.cpp  Using Constraints: Example and Test  ',' model ',
'initialize  Initialization After Constructor  ',' syntax public purpose mixed_object fixed_vec random_vec size_map example ',
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
'namespace  The CppAD::mixed Namespace Public Declarations  ',' ',
'typedef  Types Defined in the CppAD Mixed Namespace  ',' syntax begin scalar a1_double vector s_vector d_vector a1_vector sparse sparse_rc d_sparse_rcv a1_sparse_rcv end ',
'manage_gsl_rng  Set, Get, And Free A GSL Random Number Generator  ',' syntax public purpose new_gsl_rng s_in s_out get_gsl_rng free_gsl_rng example ',
'manage_gsl_rng.cpp  Manage GSL Random Number Generator: Example and Test  ',' ',
'sparse_mat_info  Sparse Matrix Information  ',' syntax public purpose row k col val resize notation sparsity pattern empty column major order lower triangular ',
'fixed_solution  Optimal Solution Returned by optimize_fixed  ',' syntax public convention fixed_opt fixed_lag bounds fix_con_lag ran_con_lag ',
'exception  CppAD Mixed Exceptions  ',' syntax public thrower brief catcher description ',
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
'whats_new_19  Changes and Additions to cppad_mixed During 2019  ',' 07-10 07-09 06-24 06-07 ',
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
