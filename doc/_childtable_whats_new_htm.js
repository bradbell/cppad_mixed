// Child table for section whats_new
document.write('\
<select onchange="whats_new_child(this)">\
<option>whats_new-&gt;</option>\
<option>whats_new_19</option>\
<option>whats_new_18</option>\
<option>whats_new_17</option>\
<option>whats_new_16</option>\
<option>whats_new_15</option>\
</select>\
');
function whats_new_child(item)
{	var child_list = [
		'whats_new_19.htm',
		'whats_new_18.htm',
		'whats_new_17.htm',
		'whats_new_16.htm',
		'whats_new_15.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
