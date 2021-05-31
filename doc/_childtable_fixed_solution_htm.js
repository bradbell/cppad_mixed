// Child table for section fixed_solution
document.write('\
<select onchange="fixed_solution_child(this)">\
<option>fixed_solution-&gt;</option>\
<option>warm_start_struct</option>\
</select>\
');
function fixed_solution_child(item)
{	var child_list = [
		'warm_start_struct.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
