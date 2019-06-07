// Child table for section fix_like_eval
document.write('\
<select onchange="fix_like_eval_child(this)">\
<option>fix_like_eval-&gt;</option>\
<option>fix_like_eval.cpp</option>\
</select>\
');
function fix_like_eval_child(item)
{	var child_list = [
		'fix_like_eval.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
