// Child table for section sparse_up_tri_sol
document.write('\
<select onchange="sparse_up_tri_sol_child(this)">\
<option>sparse_up_tri_sol-&gt;</option>\
<option>sparse_up_tri_sol.cpp</option>\
</select>\
');
function sparse_up_tri_sol_child(item)
{	var child_list = [
		'sparse_up_tri_sol.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
