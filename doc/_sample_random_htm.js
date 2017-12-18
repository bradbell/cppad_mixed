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
'sample_random.htm'
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
var list_down0 = [
'sample_random.cpp.htm'
];
var list_current0 = [
'sample_random.htm#Syntax',
'sample_random.htm#See Also',
'sample_random.htm#Prototype',
'sample_random.htm#Public',
'sample_random.htm#Purpose',
'sample_random.htm#manage_gsl_rng',
'sample_random.htm#mixed_object',
'sample_random.htm#sample',
'sample_random.htm#random_ipopt_options',
'sample_random.htm#fixed_vec',
'sample_random.htm#random_lower',
'sample_random.htm#random_upper',
'sample_random.htm#random_in',
'sample_random.htm#Covariance',
'sample_random.htm#Example'
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