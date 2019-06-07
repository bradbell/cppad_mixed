// Child table for section ran_like_hes
document.write('\
<select onchange="ran_like_hes_child(this)">\
<option>ran_like_hes-&gt;</option>\
<option>ran_like_hes.cpp</option>\
</select>\
');
function ran_like_hes_child(item)
{	var child_list = [
		'ran_like_hes.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
