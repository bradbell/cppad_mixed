// ------------------------------------------------------------ 
// Copyright (C) Bradley M. Bell 1998-2015, All rights reserved 
// ------------------------------------------------------------ 
Keyword = 
[
{ tag: 'cppad_mixed', title:'C++ Laplace Approximation of Mixed Effects Models: cppad_mixed-20210729', other:' license source code repository notation fixed theta random data z prior density p(theta) p(z|theta) p(u|theta) p(y|thetau) constraint function c(theta) optimal u^(theta) matrix a*u^(theta) problem maximum likelihood constraints negative log-density vector' },
{ tag: 'install_unix', title:'Installing cppad_mixed in Unix', other:' system requirements c++ compiler git cmake wget fortran gsl pkg-config suitesparse download special run_cmake.sh eigen ipopt paths pkg_config_path ld_library_path command check speed example installation linking using' },
{ tag: 'example_install.sh', title:'An Example Installation', other:' syntax run_test replace source' },
{ tag: 'run_cmake.sh', title:'bin/run_cmake.sh: User Configuration Options', other:' verbose_makefile build_type cmake_install_prefix debug release eigen specific_compiler extra_cxx_flags cmake_libdir cppad_mixed.pc ldlt_cholmod optimize_cppad_function for_hes_sparsity testing speed memory' },
{ tag: 'check_install.sh', title:'Example and Test Using the Installed Version of cppad_mixed', other:' syntax build_type prefixes cmake_libdir example_file pkg_config_path ld_library_path create temporary example_name main gsl_libs ipopt_libs suitesparse_libs compile link run' },
{ tag: 'theory', title:'Laplace Approximation for Mixed Effects Models', other:' reference total likelihood random f(theta u) assumption objective fixed g(theta) optimal u^(theta) h(theta r(theta) l(theta) derivative constraints approximate first order u(beta second w(beta h(beta function b(beta hessian sparse observed information' },
{ tag: 'base_class', title:'cppad_mixed: Public Declarations', other:' types user defined functions ran_likelihood fix_likelihood fix_constraint fatal_error warning constructor destructor initialize optimize_random optimize_fixed hes_fixed_obj hes_random_obj sample_fixed sample_random information_mat deprecated 2020-03-22' },
{ tag: 'derived_ctor', title:'User Defined Class Derived From cppad_mixed', other:' syntax prototype see also mixed_derived mixed_object n_fixed n_random quasi_fixed true false bool_sparsity a_rcv trace_init ... errorhandler example' },
{ tag: 'derived_ctor.cpp', title:'mixed_cppad Derived Class: Example and Test', other:'' },
{ tag: 'ran_likelihood', title:'User Defined Random Likelihood Function', other:' syntax mixed_object virtual fixed_vec random_vec constant default example' },
{ tag: 'ran_likelihood.cpp', title:'Random Likelihood: Example and Test', other:'' },
{ tag: 'fix_likelihood', title:'User Defined Fixed Likelihood Function', other:' syntax mixed_object a1_double virtual fixed_vec constant default example' },
{ tag: 'fix_likelihood.cpp', title:'Random Likelihood: Example and Test', other:'' },
{ tag: 'fix_constraint', title:'User Defined Fixed Effects Constraint Function', other:' syntax mixed_object a1_double virtual fixed_vec default example' },
{ tag: 'fix_constraint.cpp', title:'Using Constraints: Example and Test', other:' model' },
{ tag: 'initialize', title:'Initialization After Constructor', other:' syntax public purpose mixed_object fixed_vec random_vec size_map n_fixed n_random quasi_fixed a_nr a_nnz ran_like_fun.size_var fix_like_fun.size_var other fields example' },
{ tag: 'optimize_random', title:'Optimize Random Effects', other:' syntax purpose mixed_object options evaluation_method fixed_vec random_lower random_upper random_in random_out example' },
{ tag: 'optimize_random.cpp', title:'Optimize Random Effects: Example and Test', other:'' },
{ tag: 'optimize_fixed', title:'Optimize Fixed Effects', other:' syntax purpose inf mixed_object fixed_ipopt_options derivative_test hessian_approximation max_iter accept_after_max_steps nlp_scaling_method random_ipopt_options fixed_lower fixed_upper fix_constraint_lower fix_constraint_upper fixed_scale fixed_in random_lower random_upper random_in warm_start no example solution laplace ipopt_fixed' },
{ tag: 'optimize_fixed.cpp', title:'Optimize Fixed Effects: Example and Test', other:' model objective first order partials optimizer trace warm start source code' },
{ tag: 'ipopt_options', title:'An Ipopt Options Argument', other:' prototype format string integer numeric' },
{ tag: 'ipopt_trace', title:'Description of Ipopt Tracing Output', other:' discussion iter objective inf_pr inf_du lg(mu) ||d|| lg(rg) alpha_du alpha_pr ls reference' },
{ tag: 'hes_fixed_obj', title:'Compute the Hessian of The Fixed Effects Objective', other:' syntax prototype purpose mixed_object fixed_vec random_opt hes_fixed_obj_rcv example' },
{ tag: 'hes_fixed_obj.cpp', title:'Hessian of Fixed Effects Objective: Example and Test', other:' model trace_init' },
{ tag: 'hes_random_obj', title:'Compute the Hessian of The Random Effects Objective', other:' syntax prototype purpose mixed_object fixed_vec random_vec hes_random_obj_rcv example' },
{ tag: 'hes_random_obj.cpp', title:'Hessian of Random Effects Objective: Example and Test', other:'' },
{ tag: 'sample_fixed', title:'Sample Posterior for Fixed Effects', other:' syntax see also prototype purpose constant constraints covariance manage_gsl_rng mixed_object hes_fixed_obj_rcv solution fixed_lower fixed_upper error_msg example other method' },
{ tag: 'sample_fixed.cpp', title:'Sample From Fixed Effects Posterior: Example and Test', other:'' },
{ tag: 'sample_conditional', title:'Sample Posterior for Fixed Effects Using Conditional Covariance', other:' syntax replaced prototype purpose manage_gsl_rng mixed_object information_info solution fixed_lower fixed_upper random_opt theory notation subset unconstrained constraint equations example' },
{ tag: 'sample_random', title:'Simulation the Posterior Distribution for Random Effects', other:' syntax see also prototype purpose manage_gsl_rng mixed_object random_ipopt_options fixed_vec random_lower random_upper random_in covariance error_msg example' },
{ tag: 'sample_random.cpp', title:'Sample From Fixed Effects Posterior: Example and Test', other:'' },
{ tag: 'information_mat', title:'Compute the Observed Information For Fixed Effects', other:' deprecated 2020-03-22 syntax purpose mixed_object solution random_opt information_rcv example' },
{ tag: 'information_mat.cpp', title:'Observed Information Matrix: Example and Test', other:' deprecated 2020-03-22 model' },
{ tag: 'namespace', title:'The CppAD::mixed Namespace Public Declarations', other:'' },
{ tag: 'typedef', title:'Types Defined in the CppAD Mixed Namespace', other:' syntax begin scalar a1_double vector s_vector d_vector a1_vector sparse sparse_rc d_sparse_rcv a1_sparse_rcv end' },
{ tag: 'manage_gsl_rng', title:'Set, Get, And Free A GSL Random Number Generator', other:' syntax purpose new_gsl_rng s_in s_out get_gsl_rng free_gsl_rng example' },
{ tag: 'manage_gsl_rng.cpp', title:'Manage GSL Random Number Generator: Example and Test', other:'' },
{ tag: 'sparse_mat_info', title:'Sparse Matrix Information', other:' syntax purpose row k col val resize notation sparsity pattern empty column major order lower triangular' },
{ tag: 'fixed_solution', title:'Optimal Solution Returned by optimize_fixed', other:' syntax prototype convention fixed_opt fixed_lag fix_con_lag ran_con_lag warm_start trace_vec' },
{ tag: 'warm_start_struct', title:'Ipopt Warm Start Information', other:' syntax prototype scale_f x_info z_l z_u scale_x g_info lambda scale_g public' },
{ tag: 'trace_struct', title:'Ipopt Trace Information', other:' syntax prototype iter obj_value inf_pr inf_du mu d_norm regularization_size alpha_du alpha_pr ls_trials restoration' },
{ tag: 'exception', title:'CppAD Mixed Exceptions', other:' syntax thrower brief catcher description' },
{ tag: 'user_examples', title:'User API Examples', other:' functions defined called cppad_mixed programs that are speed tests demonstrate specific features' },
{ tag: 'speed', title:'Example Programs That are Also Speed and Memory Tests', other:' purpose' },
{ tag: 'ar1_xam.cpp', title:'A First Order Auto-Regressive Example and Speed Test', other:' syntax problem data p( y_t | theta ) command arguments random_seed number_random quasi_fixed trace_optimize_fixed ipopt_solve bool_sparsity hold_memory derivative_test start_near_solution output cppad_mixed_version ldlt_cholmod optimize_cppad_function ndebug_defined actual_seed initialize_bytes initialize_seconds optimize_fixed_seconds optimize_random_seconds information_mat_seconds sample_fixed_seconds final_bytes theta_0_estimate ar1_xam_ok source code' },
{ tag: 'ar1_xam.sh', title:'Example Using ar1_xam', other:' syntax test2run normal callgrind massif source code' },
{ tag: 'capture_xam.cpp', title:'A Capture Example and Speed Test', other:' syntax reference command arguments random_seed number_random quasi_fixed trace_optimize_fixed ipopt_solve bool_sparsity hold_memory derivative_test start_near_solution number_fixed_samples number_locations max_population mean_population mean_logit_probability std_logit_probability random_constraint output cppad_mixed_version ldlt_cholmod optimize_cppad_function ndebug_defined actual_seed initialize_bytes initialize_seconds optimize_fixed_seconds optimize_random_seconds information_mat_seconds sample_fixed_seconds final_bytes sum_random_effects mean_population_estimate mean_logit_probability_estimate std_logit_probability_estimate mean_population_std mean_logit_probability_std std_logit_probability_std mean_population_ratio mean_logit_probability_ratio std_probability_ratio capture_xam_ok notation p(y_it|n_iq_t) p(n_i|theta) q_t(thetau) m_i p(y_i|thetau) p(y|thetau) p(u|theta) p(theta) p(z|theta) c(theta) source code' },
{ tag: 'capture_xam.sh', title:'Example Using capture_xam', other:' syntax test2run normal callgrind massif source code' },
{ tag: 'abs_density.cpp', title:'Absolute Value In Log-Density: Example and Test', other:' model' },
{ tag: 'no_random.cpp', title:'No Random Effects: Example and Test', other:' model' },
{ tag: 'ran_constraint.cpp', title:'Constraints On Random Effects: Example and Test', other:'' },
{ tag: 'lasso.cpp', title:'Lasso on Fixed Effects: Example and Test', other:' model' },
{ tag: 'data_mismatch.cpp', title:'Random Effects Variance May Cause Data Mismatch', other:' model theory derivatives objective' },
{ tag: 'opt_ran_nan.cpp', title:'Nan\'s During Optimization of Random Effects: Example and Test', other:'' },
{ tag: 'warm_start.cpp', title:'Warm Starting Optimization: Example and Test', other:' model bounds maximum iterations optimizer trace' },
{ tag: 'release_notes', title:'Changes and Additions to cppad_mixed', other:' this year previous years' },
{ tag: 'whats_new_21', title:'Changes and Additions to cppad_mixed During 2021', other:' 07-29 07-27 07-13 06-29 06-26 06-24 06-13 06-10 06-08 06-07 06-06 06-03 05-31 05-29 05-28 05-27 05-26 05-17 05-15 05-11 05-07 05-01 03-02' },
{ tag: 'whats_new_20', title:'Changes and Additions to cppad_mixed During 2020', other:' 12-22 11-30 11-23 11-21 11-19 11-18 11-05 11-04 11-03 11-02 10-31 10-21 10-10 10-06 08-31 08-21 07-02 06-30 06-07 05-30 05-29 05-27 03-28 03-25 03-23 03-22 03-18 03-15 api' },
{ tag: 'whats_new_19', title:'Changes and Additions to cppad_mixed During 2019', other:' 10-08 09-30 07-24 07-23 07-20 07-19 07-10 07-09 06-24 06-07' },
{ tag: 'whats_new_18', title:'Changes and Additions to cppad_mixed During 2018', other:' 10-08 09-25 08-30 08-27 08-22 08-20 08-18 08-09 08-08 07-25 07-12 06-30 06-29 06-20 06-14 06-09 06-04 05-21 05-07 05-03 04-06 03-22 03-10 02-20 02-12 02-11 02-10 02-08 02-07 02-05 01-23 01-22 01-21 01-15 01-14' },
{ tag: 'whats_new_17', title:'Changes and Additions to cppad_mixed During 2017', other:' 12-28 12-25 12-22 12-16 12-10 10-27 10-24 10-09 10-07 09-30 09-23 09-21 09-18 09-16 09-15 09-14 09-02 08-30 08-01 04-24 04-23 04-06 04-02 03-27 03-25 03-23 03-20 03-12 03-11 03-10 03-09 03-08 03-06 03-02 03-01 01-26 01-24 01-22 01-14' },
{ tag: 'whats_new_16', title:'Changes and Additions to cppad_mixed During 2016', other:' 11-09 11-07 11-04 11-02 10-30 10-28 10-27 10-25 10-24 10-18 10-16 10-15 10-14 10-06 09-30 09-27 09-24 09-23 07-28 07-27 07-26 07-25 07-20 07-18 07-14 07-13 07-12 07-10 07-09 06-24 06-22 06-19 06-18 06-17 06-13 06-12 06-11 06-07 06-06 06-05 06-04 06-03 optimize_random optimize_fixed sample_random 05-15 05-11 05-08 05-06 05-04 05-03 04-29 04-27 04-23 04-19 04-18 04-17 04-16 04-10 04-15 04-09 04-08 04-07 04-06 04-05 04-03 04-02 04-01 03-29 03-28 03-09 02-26 02-06 01-26 01-25 01-22 01-21 01-19 01-16 01-15 01-14 01-13 01-10 01-09 01-05 01-04 01-01' },
{ tag: 'whats_new_15', title:'Changes and Additions to cppad_mixed During 2015', other:' 12-25 12-24 12-16 12-14 12-13 12-10' },
{ tag: 'wish_list', title:'CppAD Mixed Wish List', other:' better asymptotic statistics multi-threading fixed likelihood hessian examples d_sparse_rcv aborting optimization windows install' },
{ tag: 'math_notation', title:'Mathematical Notation', other:' b c_l c_u f g p r u^(theta) w y z' }
]

var MaxList = 100;
var Nstring = -1;
var Nkeyword = Keyword.length;
var Row2Tag  = [];
Initialize();

function Initialize()
{
	UpdateList();
	document.search.keywords.focus();
}
function UpdateList()
{
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

	var nlist       = 0;
	var title_list  = '';
	var tag_list    = '';
	for(i = 0; (i < Nkeyword) && (nlist < MaxList) ; i++)
	{
		var match = true;
		for(j = 0; j < nword; j++)
		{	var flag = pattern[j].test(Keyword[i].tag);
			flag     = flag || pattern[j].test(Keyword[i].title);
			flag     = flag || pattern[j].test(Keyword[i].other);
			match    = match && flag;
		}
		if( match )
		{
			var tag    = Keyword[i].tag;
			var title  = Keyword[i].title
			title      = title.split(/\s+/);
			title      = title.join(' ');
			title_list = title_list + title + '\n';
			tag_list   = tag_list + tag + '\n'
			Row2Tag[nlist] = tag;
			nlist = nlist + 1;
		}
	}
	document.search.title_list.value = title_list;
	document.search.title_list.setAttribute('wrap', 'off');;
	document.search.title_list.readOnly = true;
	document.search.tag_list.value = tag_list;
	document.search.tag_list.setAttribute('wrap', 'off');;
	document.search.tag_list.readOnly = true;
}
function Choose(textarea)
{	var start_select = textarea.value.substring(0, textarea.selectionStart);
	var start_pos    = Math.max(0, start_select.lastIndexOf('\n') + 1);
	var length       = textarea.value.length;
	var end_select   = 
		textarea.value.substring(textarea.selectionEnd, length);
	var end_pos  = end_select.indexOf('\n');
	end_pos      = textarea.selectionEnd + end_pos;
	textarea.selectionStart = start_pos;
	textarea.selectionEnd   = end_pos;
	var row = start_select.split('\n').length - 1;
	var tag = Row2Tag[row];
	document.search.selection.value    = tag.toLowerCase();
	document.search.selection.readOnly = true;
	
	return true;
}
function Goto()
{  selection       = document.search.selection.value;
   if( selection != '' )
	    parent.location = selection + '.htm';
}
function CheckForReturn()
{
	var key = event.which;
	if( key == 13 ) Goto();
}