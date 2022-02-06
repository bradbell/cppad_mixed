// Child table for section release_notes
document.write('\
<select onchange="release_notes_child(this)">\
<option>release_notes-&gt;</option>\
<option>whats_new_22</option>\
<option>whats_new_21</option>\
<option>whats_new_20</option>\
<option>whats_new_19</option>\
<option>whats_new_18</option>\
<option>whats_new_17</option>\
<option>whats_new_16</option>\
<option>whats_new_15</option>\
</select>\
');
function release_notes_child(item)
{	var child_list = [
		'whats_new_22.htm',
		'whats_new_21.htm',
		'whats_new_20.htm',
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
