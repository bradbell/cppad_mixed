// Child table for section sparse_mat2low
document.write('\
<select onchange="sparse_mat2low_child(this)">\
<option>sparse_mat2low-&gt;</option>\
<option>sparse_mat2low.cpp</option>\
</select>\
');
function sparse_mat2low_child(item)
{	var child_list = [
		'sparse_mat2low.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
