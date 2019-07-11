// Child table for section namespace
document.write('\
<select onchange="namespace_child(this)">\
<option>namespace-&gt;</option>\
<option>typedef</option>\
<option>manage_gsl_rng</option>\
<option>sparse_mat_info</option>\
<option>fixed_solution</option>\
<option>exception</option>\
</select>\
');
function namespace_child(item)
{	var child_list = [
		'typedef.htm',
		'manage_gsl_rng.htm',
		'sparse_mat_info.htm',
		'fixed_solution.htm',
		'exception.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
