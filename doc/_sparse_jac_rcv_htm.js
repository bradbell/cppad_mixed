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
'sparse_jac_rcv.htm'
];
var list_down2 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_17.htm',
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
var list_current0 = [
'sparse_jac_rcv.htm#Syntax',
'sparse_jac_rcv.htm#Private',
'sparse_jac_rcv.htm#Purpose',
'sparse_jac_rcv.htm#subset',
'sparse_jac_rcv.htm#subset.nnz',
'sparse_jac_rcv.htm#subset.row',
'sparse_jac_rcv.htm#subset.col',
'sparse_jac_rcv.htm#subset.val',
'sparse_jac_rcv.htm#forward',
'sparse_jac_rcv.htm#work',
'sparse_jac_rcv.htm#Computing Sparse Jacobians',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.Forward Mode',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.Reverse Mode',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.f',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.group_max',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.x',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.not_used_pattern',
'sparse_jac_rcv.htm#Computing Sparse Jacobians.not_used_coloring'
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
