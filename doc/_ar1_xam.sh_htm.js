var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'user.htm',
'speed.htm',
'ar1_xam.cpp.htm',
'ar1_xam.sh.htm'
];
var list_down4 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_17.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down3 = [
'speed.htm',
'abs_density.cpp.htm',
'no_random.cpp.htm',
'ran_constraint.cpp.htm',
'lasso.cpp.htm',
'data_mismatch.cpp.htm',
'opt_ran_nan.cpp.htm'
];
var list_down2 = [
'ar1_xam.cpp.htm',
'capture_xam.cpp.htm'
];
var list_down1 = [
'ar1_xam.sh.htm'
];
var list_current0 = [
'ar1_xam.sh.htm#Syntax',
'ar1_xam.sh.htm#test2run',
'ar1_xam.sh.htm#test2run.normal',
'ar1_xam.sh.htm#test2run.callgrind',
'ar1_xam.sh.htm#test2run.massif',
'ar1_xam.sh.htm#test2run.Source Code'
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
function choose_down4(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down4[index-1];
}
function choose_down3(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down3[index-1];
}
function choose_down2(item)
{	var index          = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = list_down2[index-1];
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
