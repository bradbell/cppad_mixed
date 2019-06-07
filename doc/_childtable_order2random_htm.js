// Child table for section order2random
document.write('\
<select onchange="order2random_child(this)">\
<option>order2random-&gt;</option>\
<option>order2random.cpp</option>\
</select>\
');
function order2random_child(item)
{	var child_list = [
		'order2random.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
