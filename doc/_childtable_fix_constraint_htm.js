// Child table for section fix_constraint
document.write('\
<select onchange="fix_constraint_child(this)">\
<option>fix_constraint-&gt;</option>\
<option>fix_constraint.cpp</option>\
</select>\
');
function fix_constraint_child(item)
{	var child_list = [
		'fix_constraint.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
