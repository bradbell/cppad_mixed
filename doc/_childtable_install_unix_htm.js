// Child table for section install_unix
document.write('\
<select onchange="install_unix_child(this)">\
<option>install_unix-&gt;</option>\
<option>example_install.sh</option>\
<option>run_cmake.sh</option>\
<option>check_install.sh</option>\
</select>\
');
function install_unix_child(item)
{	var child_list = [
		'example_install.sh.htm',
		'run_cmake.sh.htm',
		'check_install.sh.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
