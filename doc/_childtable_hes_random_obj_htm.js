// Child table for section hes_random_obj
document.write('\
<select onchange="hes_random_obj_child(this)">\
<option>hes_random_obj-&gt;</option>\
<option>hes_random_obj.cpp</option>\
</select>\
');
function hes_random_obj_child(item)
{	var child_list = [
		'hes_random_obj.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
