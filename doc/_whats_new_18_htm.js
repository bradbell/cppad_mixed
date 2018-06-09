var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'whats_new_18.htm'
];
var list_down1 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_18.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down0 = [
'whats_new_17.htm',
'whats_new_16.htm',
'whats_new_15.htm'
];
var list_current0 = [
'whats_new_18.htm#Contents',
'whats_new_18.htm#06-09',
'whats_new_18.htm#06-04',
'whats_new_18.htm#05-21',
'whats_new_18.htm#05-07',
'whats_new_18.htm#05-03',
'whats_new_18.htm#04-06',
'whats_new_18.htm#03-22',
'whats_new_18.htm#03-10',
'whats_new_18.htm#02-20',
'whats_new_18.htm#02-12',
'whats_new_18.htm#02-11',
'whats_new_18.htm#02-10',
'whats_new_18.htm#02-08',
'whats_new_18.htm#02-07',
'whats_new_18.htm#02-05',
'whats_new_18.htm#01-23',
'whats_new_18.htm#01-22',
'whats_new_18.htm#01-21',
'whats_new_18.htm#01-15',
'whats_new_18.htm#01-14'
];
function choose_across0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_across0[index-1];
}
function choose_up0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_up0[index-1];
}
function choose_down1(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down1[index-1];
}
function choose_down0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down0[index-1];
}
function choose_current0(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_current0[index-1];
}
