// Child table for section optimize_random
document.write('\
<select onchange="optimize_random_child(this)">\
<option>optimize_random-&gt;</option>\
<option>optimize_random.cpp</option>\
</select>\
');
function optimize_random_child(item)
{	var child_list = [
		'optimize_random.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
