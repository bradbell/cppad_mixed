// Child table for section public
document.write('\
<select onchange="public_child(this)">\
<option>public-&gt;</option>\
<option>derived_ctor</option>\
<option>ran_likelihood</option>\
<option>fix_likelihood</option>\
<option>fix_constraint</option>\
<option>initialize</option>\
<option>optimize_random</option>\
<option>optimize_fixed</option>\
<option>information_mat</option>\
<option>sample_fixed</option>\
<option>sample_random</option>\
</select>\
');
function public_child(item)
{	var child_list = [
		'derived_ctor.htm',
		'ran_likelihood.htm',
		'fix_likelihood.htm',
		'fix_constraint.htm',
		'initialize.htm',
		'optimize_random.htm',
		'optimize_fixed.htm',
		'information_mat.htm',
		'sample_fixed.htm',
		'sample_random.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
