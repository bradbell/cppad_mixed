// Child table for section ran_obj_eval
document.write('\
<select onchange="ran_obj_eval_child(this)">\
<option>ran_obj_eval-&gt;</option>\
<option>ran_obj_eval.cpp</option>\
</select>\
');
function ran_obj_eval_child(item)
{	var child_list = [
		'ran_obj_eval.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
