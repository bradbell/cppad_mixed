// Child table for section sparse_rcv2eigen
document.write('\
<select onchange="sparse_rcv2eigen_child(this)">\
<option>sparse_rcv2eigen-&gt;</option>\
<option>sparse_rcv2eigen.cpp</option>\
</select>\
');
function sparse_rcv2eigen_child(item)
{	var child_list = [
		'sparse_rcv2eigen.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
