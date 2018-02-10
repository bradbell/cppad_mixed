var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'namespace.htm',
'newton_step.htm'
];
var list_down2 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_18.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down1 = [
'typedef.htm',
'configure.hpp.htm',
'exception.htm',
'fixed_solution.htm',
'ipopt_app_status.htm',
'ipopt_fixed.htm',
'ipopt_random.htm',
'ldlt_cholmod.htm',
'ldlt_eigen.htm',
'manage_gsl_rng.htm',
'newton_step.htm',
'sparse_hes_rcv.htm',
'sparse_hes_info.htm',
'sparse_jac_rcv.htm',
'sparse_mat_info.htm',
'triple2eigen.htm',
'undetermined.htm',
'sparse_low_tri_sol.htm',
'sparse_up_tri_sol.htm',
'sparse_scale_diag.htm',
'sparse_low2sym.htm',
'sparse_mat2low.htm',
'sparse_eigen2info.htm',
'sparse_info2eigen.htm',
'sparse_print.htm',
'sparsity_print.htm',
'sparse_ad_cholesky.htm'
];
var list_down0 = [
'newton_step.cpp.htm',
'newton_step_algo_ctor.htm',
'newton_step_algo.htm',
'newton_step_ctor.htm',
'newton_step_initialize.htm',
'newton_step_size_var.htm',
'newton_step_eval.htm'
];
var list_current0 = [
'newton_step.htm#Syntax',
'newton_step.htm#Private',
'newton_step.htm#Purpose',
'newton_step.htm#Constructor',
'newton_step.htm#Destructor',
'newton_step.htm#initialize',
'newton_step.htm#initialize.a1fun',
'newton_step.htm#initialize.hes_rcv',
'newton_step.htm#theta',
'newton_step.htm#u',
'newton_step.htm#size_var',
'newton_step.htm#eval',
'newton_step.htm#eval.a1_theta_u_v',
'newton_step.htm#eval.a1_logdet_step',
'newton_step.htm#eval.Checkpoint',
'newton_step.htm#Example',
'newton_step.htm#Contents'
];
function choose_across0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_across0[index-1];
}
function choose_up0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_up0[index-1];
}
function choose_down2(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down2[index-1];
}
function choose_down1(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down1[index-1];
}
function choose_down0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down0[index-1];
}
function choose_current0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_current0[index-1];
}
