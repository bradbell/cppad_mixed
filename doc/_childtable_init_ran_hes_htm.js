// Child table for section init_ran_hes
document.write('\
<select onchange="init_ran_hes_child(this)">\
<option>init_ran_hes-&gt;</option>\
<option>ran_hes_fun.cpp</option>\
</select>\
');
function init_ran_hes_child(item)
{	var child_list = [
		'ran_hes_fun.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
