// Child table for section undetermined
document.write('\
<select onchange="undetermined_child(this)">\
<option>undetermined-&gt;</option>\
<option>undetermined.cpp</option>\
</select>\
');
function undetermined_child(item)
{	var child_list = [
		'undetermined.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
