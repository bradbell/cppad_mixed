// Child table for section sparse_ad_cholesky
document.write('\
<select onchange="sparse_ad_cholesky_child(this)">\
<option>sparse_ad_cholesky-&gt;</option>\
<option>sparse_ad_cholesky_initialize</option>\
<option>sparse_ad_cholesky_p</option>\
<option>sparse_ad_cholesky_eval</option>\
<option>set_jac_sparsity</option>\
<option>set_hes_sparsity</option>\
<option>sparse_ad_chol_eval.cpp</option>\
<option>sparse_ad_chol_perm.cpp</option>\
<option>sparse_ad_chol_eq.cpp</option>\
<option>sparse_ad_chol_var.cpp</option>\
<option>sparse_ad_chol_sp1.cpp</option>\
<option>sparse_ad_chol_sp2.cpp</option>\
</select>\
');
function sparse_ad_cholesky_child(item)
{	var child_list = [
		'sparse_ad_cholesky_initialize.htm',
		'sparse_ad_cholesky_p.htm',
		'sparse_ad_cholesky_eval.htm',
		'set_jac_sparsity.htm',
		'set_hes_sparsity.htm',
		'sparse_ad_chol_eval.cpp.htm',
		'sparse_ad_chol_perm.cpp.htm',
		'sparse_ad_chol_eq.cpp.htm',
		'sparse_ad_chol_var.cpp.htm',
		'sparse_ad_chol_sp1.cpp.htm',
		'sparse_ad_chol_sp2.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
