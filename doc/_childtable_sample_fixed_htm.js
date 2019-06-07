// Child table for section sample_fixed
document.write('\
<select onchange="sample_fixed_child(this)">\
<option>sample_fixed-&gt;</option>\
<option>sample_fixed.cpp</option>\
<option>sample_conditional</option>\
</select>\
');
function sample_fixed_child(item)
{	var child_list = [
		'sample_fixed.cpp.htm',
		'sample_conditional.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
