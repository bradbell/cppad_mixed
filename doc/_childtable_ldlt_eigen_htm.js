// Child table for section ldlt_eigen
document.write('\
<select onchange="ldlt_eigen_child(this)">\
<option>ldlt_eigen-&gt;</option>\
<option>ldlt_eigen_ctor</option>\
<option>ldlt_eigen_init</option>\
<option>ldlt_eigen_pattern</option>\
<option>ldlt_eigen_update</option>\
<option>ldlt_eigen_split</option>\
<option>ldlt_eigen_logdet</option>\
<option>ldlt_eigen_solve_H</option>\
<option>ldlt_eigen_sim_cov</option>\
<option>ldlt_eigen_inv</option>\
<option>ldlt_eigen_solve_LDLT</option>\
<option>ldlt_eigen.cpp</option>\
</select>\
');
function ldlt_eigen_child(item)
{	var child_list = [
		'ldlt_eigen_ctor.htm',
		'ldlt_eigen_init.htm',
		'ldlt_eigen_pattern.htm',
		'ldlt_eigen_update.htm',
		'ldlt_eigen_split.htm',
		'ldlt_eigen_logdet.htm',
		'ldlt_eigen_solve_h.htm',
		'ldlt_eigen_sim_cov.htm',
		'ldlt_eigen_inv.htm',
		'ldlt_eigen_solve_ldlt.htm',
		'ldlt_eigen.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
