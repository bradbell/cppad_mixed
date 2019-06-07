// Child table for section init_ran_jac
document.write('\
<select onchange="init_ran_jac_child(this)">\
<option>init_ran_jac-&gt;</option>\
<option>ran_jac_fun.cpp</option>\
</select>\
');
function init_ran_jac_child(item)
{	var child_list = [
		'ran_jac_fun.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
