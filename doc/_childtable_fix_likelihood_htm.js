// Child table for section fix_likelihood
document.write('\
<select onchange="fix_likelihood_child(this)">\
<option>fix_likelihood-&gt;</option>\
<option>fix_likelihood.cpp</option>\
</select>\
');
function fix_likelihood_child(item)
{	var child_list = [
		'fix_likelihood.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
