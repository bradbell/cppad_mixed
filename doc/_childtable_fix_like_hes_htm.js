// Child table for section fix_like_hes
document.write('\
<select onchange="fix_like_hes_child(this)">\
<option>fix_like_hes-&gt;</option>\
<option>fix_like_hes.cpp</option>\
</select>\
');
function fix_like_hes_child(item)
{	var child_list = [
		'fix_like_hes.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
