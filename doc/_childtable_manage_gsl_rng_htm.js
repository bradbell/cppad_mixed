// Child table for section manage_gsl_rng
document.write('\
<select onchange="manage_gsl_rng_child(this)">\
<option>manage_gsl_rng-&gt;</option>\
<option>manage_gsl_rng.cpp</option>\
</select>\
');
function manage_gsl_rng_child(item)
{	var child_list = [
		'manage_gsl_rng.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
