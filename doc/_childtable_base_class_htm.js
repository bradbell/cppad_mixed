// Child table for section base_class
document.write('\
<select onchange="base_class_child(this)">\
<option>base_class-&gt;</option>\
<option>public</option>\
<option>private</option>\
</select>\
');
function base_class_child(item)
{	var child_list = [
		'public.htm',
		'private.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
