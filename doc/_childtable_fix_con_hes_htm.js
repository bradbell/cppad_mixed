// Child table for section fix_con_hes
document.write('\
<select onchange="fix_con_hes_child(this)">\
<option>fix_con_hes-&gt;</option>\
<option>fix_con_hes.cpp</option>\
</select>\
');
function fix_con_hes_child(item)
{	var child_list = [
		'fix_con_hes.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
