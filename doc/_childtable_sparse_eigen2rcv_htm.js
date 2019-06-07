// Child table for section sparse_eigen2rcv
document.write('\
<select onchange="sparse_eigen2rcv_child(this)">\
<option>sparse_eigen2rcv-&gt;</option>\
<option>sparse_eigen2rcv.cpp</option>\
</select>\
');
function sparse_eigen2rcv_child(item)
{	var child_list = [
		'sparse_eigen2rcv.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
