// Child table for section speed
document.write('\
<select onchange="speed_child(this)">\
<option>speed-&gt;</option>\
<option>ar1_xam.cpp</option>\
<option>capture_xam.cpp</option>\
</select>\
');
function speed_child(item)
{	var child_list = [
		'ar1_xam.cpp.htm',
		'capture_xam.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
