// Child table for section namespace
document.write('\
<select onchange="namespace_child(this)">\
<option>namespace-&gt;</option>\
<option>typedef</option>\
<option>configure.hpp</option>\
<option>exception</option>\
<option>fixed_solution</option>\
<option>ipopt_app_status</option>\
<option>ipopt_fixed</option>\
<option>ipopt_random</option>\
<option>ldlt_cholmod</option>\
<option>ldlt_eigen</option>\
<option>manage_gsl_rng</option>\
<option>sparse_hes_rcv</option>\
<option>sparse_hes_info</option>\
<option>sparse_jac_rcv</option>\
<option>sparse_mat_info</option>\
<option>triple2eigen</option>\
<option>order2random</option>\
<option>undetermined</option>\
<option>sparse_low_tri_sol</option>\
<option>sparse_up_tri_sol</option>\
<option>sparse_scale_diag</option>\
<option>sparse_low2sym</option>\
<option>sparse_mat2low</option>\
<option>sparse_eigen2info</option>\
<option>sparse_info2eigen</option>\
<option>sparse_eigen2rcv</option>\
<option>sparse_rcv2eigen</option>\
<option>sparse_print</option>\
<option>sparsity_print</option>\
<option>sparse_ad_cholesky</option>\
</select>\
');
function namespace_child(item)
{	var child_list = [
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
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
