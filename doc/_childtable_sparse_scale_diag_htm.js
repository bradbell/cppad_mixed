// Child table for section sparse_scale_diag
document.write('\
<select onchange="sparse_scale_diag_child(this)">\
<option>sparse_scale_diag-&gt;</option>\
<option>sparse_scale_diag.cpp</option>\
</select>\
');
function sparse_scale_diag_child(item)
{	var child_list = [
		'sparse_scale_diag.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
