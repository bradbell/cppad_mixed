var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'_reference.htm'
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
var list_current0 = [
'_reference.htm#A',
'_reference.htm#B',
'_reference.htm#C',
'_reference.htm#D',
'_reference.htm#E',
'_reference.htm#F',
'_reference.htm#H',
'_reference.htm#I',
'_reference.htm#L',
'_reference.htm#M',
'_reference.htm#N',
'_reference.htm#O',
'_reference.htm#P',
'_reference.htm#R',
'_reference.htm#S',
'_reference.htm#T',
'_reference.htm#U',
'_reference.htm#W'
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
