// Child table for section update_factor
document.write('\
<select onchange="update_factor_child(this)">\
<option>update_factor-&gt;</option>\
<option>update_factor.cpp</option>\
</select>\
');
function update_factor_child(item)
{	var child_list = [
		'update_factor.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
