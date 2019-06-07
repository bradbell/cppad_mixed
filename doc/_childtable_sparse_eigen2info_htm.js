// Child table for section sparse_eigen2info
document.write('\
<select onchange="sparse_eigen2info_child(this)">\
<option>sparse_eigen2info-&gt;</option>\
<option>sparse_eigen2info.cpp</option>\
</select>\
');
function sparse_eigen2info_child(item)
{	var child_list = [
		'sparse_eigen2info.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
