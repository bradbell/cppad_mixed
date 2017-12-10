var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'base_class.htm',
'public.htm',
'initialize.htm'
];
var list_down3 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_17.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down2 = [
'public.htm',
'private.htm'
];
var list_down1 = [
'derived_ctor.htm',
'ran_likelihood.htm',
'fix_likelihood.htm',
'fix_constraint.htm',
'initialize.htm',
'optimize_random.htm',
'optimize_fixed.htm',
'information_mat.htm',
'sample_fixed.htm',
'sample_random.htm',
'ran_likelihood_jac.htm',
'ran_likelihood_hes.htm'
];
var list_current0 = [
'initialize.htm#Syntax',
'initialize.htm#Public',
'initialize.htm#Purpose',
'initialize.htm#mixed_object',
'initialize.htm#fixed_vec',
'initialize.htm#random_vec',
'initialize.htm#ran_likelihood_jac',
'initialize.htm#ran_likelihood_hes',
'initialize.htm#size_map',
'initialize.htm#Example'
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
