// Child table for section sample_random
document.write('\
<select onchange="sample_random_child(this)">\
<option>sample_random-&gt;</option>\
<option>sample_random.cpp</option>\
</select>\
');
function sample_random_child(item)
{	var child_list = [
		'sample_random.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
