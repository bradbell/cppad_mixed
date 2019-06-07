// Child table for section optimize_fixed
document.write('\
<select onchange="optimize_fixed_child(this)">\
<option>optimize_fixed-&gt;</option>\
<option>optimize_fixed.cpp</option>\
<option>ipopt_options</option>\
<option>ipopt_trace</option>\
</select>\
');
function optimize_fixed_child(item)
{	var child_list = [
		'optimize_fixed.cpp.htm',
		'ipopt_options.htm',
		'ipopt_trace.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
