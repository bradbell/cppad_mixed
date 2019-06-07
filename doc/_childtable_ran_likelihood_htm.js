// Child table for section ran_likelihood
document.write('\
<select onchange="ran_likelihood_child(this)">\
<option>ran_likelihood-&gt;</option>\
<option>ran_likelihood.cpp</option>\
</select>\
');
function ran_likelihood_child(item)
{	var child_list = [
		'ran_likelihood.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
