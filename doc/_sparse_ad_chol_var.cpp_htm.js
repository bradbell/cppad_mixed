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
'sparse_ad_cholesky.htm',
'sparse_ad_chol_var.cpp.htm'
];
var list_down3 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_17.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down2 = [
'typedef.htm',
'configure.hpp.htm',
'exception.htm',
'fixed_solution.htm',
'ipopt_app_status.htm',
'ipopt_fixed.htm',
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
'sparse_ad_cholesky.htm',
'sparse_print.htm',
'sparsity_print.htm'
];
var list_down1 = [
'sparse_ad_cholesky_initialize.htm',
'sparse_ad_cholesky_p.htm',
'sparse_ad_cholesky_eval.htm',
'set_jac_sparsity.htm',
'set_hes_sparsity.htm',
'sparse_ad_chol_eval.cpp.htm',
'sparse_ad_chol_perm.cpp.htm',
'sparse_ad_chol_eq.cpp.htm',
'sparse_ad_chol_var.cpp.htm',
'sparse_ad_chol_sp.cpp.htm'
];
var list_current0 = [
'sparse_ad_chol_var.cpp.htm#Problem',
'sparse_ad_chol_var.cpp.htm#Permutation',
'sparse_ad_chol_var.cpp.htm#Cholesky Factor',
'sparse_ad_chol_var.cpp.htm#Source'
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
function choose_down3(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down3[index-1];
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