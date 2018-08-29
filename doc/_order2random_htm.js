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
'order2random.htm'
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
'sparse_hes_rcv.htm',
'sparse_hes_info.htm',
'sparse_jac_rcv.htm',
'sparse_mat_info.htm',
'triple2eigen.htm',
'order2random.htm',
'undetermined.htm',
'sparse_low_tri_sol.htm',
'sparse_up_tri_sol.htm',
'sparse_scale_diag.htm',
'sparse_low2sym.htm',
'sparse_mat2low.htm',
'sparse_eigen2info.htm',
'sparse_info2eigen.htm',
'sparse_eigen2rcv.htm',
'sparse_rcv2eigen.htm',
'sparse_print.htm',
'sparsity_print.htm',
'sparse_ad_cholesky.htm'
];
var list_down0 = [
'order2random.cpp.htm'
];
var list_current0 = [
'order2random.htm#Syntax',
'order2random.htm#Prototype',
'order2random.htm#Private',
'order2random.htm#mixed_object',
'order2random.htm#n_fixed',
'order2random.htm#n_random',
'order2random.htm#jac_a1fun',
'order2random.htm#a1_ldlt_ran_hes',
'order2random.htm#beta',
'order2random.htm#theta_u',
'order2random.htm#W',
'order2random.htm#Example'
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
