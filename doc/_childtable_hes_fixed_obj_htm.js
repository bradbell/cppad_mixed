// Child table for section hes_fixed_obj
document.write('\
<select onchange="hes_fixed_obj_child(this)">\
<option>hes_fixed_obj-&gt;</option>\
<option>hes_fixed_obj.cpp</option>\
</select>\
');
function hes_fixed_obj_child(item)
{	var child_list = [
		'hes_fixed_obj.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
