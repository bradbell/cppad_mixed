// Child table for section cppad_mixed
document.write('\
<select onchange="cppad_mixed_child(this)">\
<option>cppad_mixed-&gt;</option>\
<option>install_unix</option>\
<option>theory</option>\
<option>public</option>\
<option>public_namespace</option>\
<option>user</option>\
<option>whats_new</option>\
<option>wish_list</option>\
<option>math_notation</option>\
<option>_reference</option>\
<option>_index</option>\
<option>_search</option>\
<option>_external</option>\
</select>\
');
function cppad_mixed_child(item)
{	var child_list = [
		'install_unix.htm',
		'theory.htm',
		'public.htm',
		'public_namespace.htm',
		'user.htm',
		'whats_new.htm',
		'wish_list.htm',
		'math_notation.htm',
		'_reference.htm',
		'_index.htm',
		'_search.htm',
		'_external.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
