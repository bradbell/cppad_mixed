// Child table for section laplace_obj_hes
document.write('\
<select onchange="laplace_obj_hes_child(this)">\
<option>laplace_obj_hes-&gt;</option>\
<option>laplace_obj_hes.cpp</option>\
</select>\
');
function laplace_obj_hes_child(item)
{	var child_list = [
		'laplace_obj_hes.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
