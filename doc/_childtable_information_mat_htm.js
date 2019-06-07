// Child table for section information_mat
document.write('\
<select onchange="information_mat_child(this)">\
<option>information_mat-&gt;</option>\
<option>information_mat.cpp</option>\
</select>\
');
function information_mat_child(item)
{	var child_list = [
		'information_mat.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
