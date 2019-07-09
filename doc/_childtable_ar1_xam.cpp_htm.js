// Child table for section ar1_xam.cpp
document.write('\
<select onchange="ar1_xam__46__cpp_child(this)">\
<option>ar1_xam.cpp-&gt;</option>\
<option>ar1_xam.sh</option>\
</select>\
');
function ar1_xam__46__cpp_child(item)
{	var child_list = [
		'ar1_xam.sh.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
