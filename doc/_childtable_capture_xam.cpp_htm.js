// Child table for section capture_xam.cpp
document.write('\
<select onchange="capture_xam__46__cpp_child(this)">\
<option>capture_xam.cpp-&gt;</option>\
<option>capture_xam.sh</option>\
</select>\
');
function capture_xam__46__cpp_child(item)
{	var child_list = [
		'capture_xam.sh.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
