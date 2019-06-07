// Child table for section derived_ctor
document.write('\
<select onchange="derived_ctor_child(this)">\
<option>derived_ctor-&gt;</option>\
<option>derived_ctor.cpp</option>\
</select>\
');
function derived_ctor_child(item)
{	var child_list = [
		'derived_ctor.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
