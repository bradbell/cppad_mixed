// Child table for section init_hes_cross
document.write('\
<select onchange="init_hes_cross_child(this)">\
<option>init_hes_cross-&gt;</option>\
<option>hes_cross.cpp</option>\
</select>\
');
function init_hes_cross_child(item)
{	var child_list = [
		'hes_cross.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
