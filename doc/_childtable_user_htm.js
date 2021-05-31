// Child table for section user
document.write('\
<select onchange="user_child(this)">\
<option>user-&gt;</option>\
<option>speed</option>\
<option>abs_density.cpp</option>\
<option>no_random.cpp</option>\
<option>ran_constraint.cpp</option>\
<option>lasso.cpp</option>\
<option>data_mismatch.cpp</option>\
<option>opt_ran_nan.cpp</option>\
<option>warm_start.cpp</option>\
</select>\
');
function user_child(item)
{	var child_list = [
		'speed.htm',
		'abs_density.cpp.htm',
		'no_random.cpp.htm',
		'ran_constraint.cpp.htm',
		'lasso.cpp.htm',
		'data_mismatch.cpp.htm',
		'opt_ran_nan.cpp.htm',
		'warm_start.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
