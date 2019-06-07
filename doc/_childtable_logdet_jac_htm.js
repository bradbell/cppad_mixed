// Child table for section logdet_jac
document.write('\
<select onchange="logdet_jac_child(this)">\
<option>logdet_jac-&gt;</option>\
<option>logdet_jac.cpp</option>\
</select>\
');
function logdet_jac_child(item)
{	var child_list = [
		'logdet_jac.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
