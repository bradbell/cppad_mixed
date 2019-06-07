// Child table for section sparse_low2sym
document.write('\
<select onchange="sparse_low2sym_child(this)">\
<option>sparse_low2sym-&gt;</option>\
<option>sparse_low2sym.cpp</option>\
</select>\
');
function sparse_low2sym_child(item)
{	var child_list = [
		'sparse_low2sym.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
