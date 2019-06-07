// Child table for section ran_obj_jac
document.write('\
<select onchange="ran_obj_jac_child(this)">\
<option>ran_obj_jac-&gt;</option>\
<option>ran_obj_jac.cpp</option>\
</select>\
');
function ran_obj_jac_child(item)
{	var child_list = [
		'ran_obj_jac.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
