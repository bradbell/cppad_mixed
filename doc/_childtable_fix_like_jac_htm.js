// Child table for section fix_like_jac
document.write('\
<select onchange="fix_like_jac_child(this)">\
<option>fix_like_jac-&gt;</option>\
<option>fix_like_jac.cpp</option>\
</select>\
');
function fix_like_jac_child(item)
{	var child_list = [
		'fix_like_jac.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
