var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'math_notation.htm'
];
var list_down1 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_17.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_current0 = [
'math_notation.htm#A',
'math_notation.htm#B',
'math_notation.htm#c',
'math_notation.htm#c_L',
'math_notation.htm#c_U',
'math_notation.htm#f',
'math_notation.htm#g',
'math_notation.htm#h',
'math_notation.htm#H',
'math_notation.htm#p',
'math_notation.htm#r',
'math_notation.htm#L',
'math_notation.htm#u',
'math_notation.htm#u^(theta)',
'math_notation.htm#U',
'math_notation.htm#W',
'math_notation.htm#theta',
'math_notation.htm#y',
'math_notation.htm#z'
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
