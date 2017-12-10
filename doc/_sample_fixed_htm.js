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
'sample_fixed.htm'
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
'sample_fixed.cpp.htm',
'sample_conditional.htm'
];
var list_current0 = [
'sample_fixed.htm#Syntax',
'sample_fixed.htm#See Also',
'sample_fixed.htm#Prototype',
'sample_fixed.htm#Public',
'sample_fixed.htm#Purpose',
'sample_fixed.htm#manage_gsl_rng',
'sample_fixed.htm#mixed_object',
'sample_fixed.htm#sample',
'sample_fixed.htm#information_rcv',
'sample_fixed.htm#solution',
'sample_fixed.htm#fixed_lower',
'sample_fixed.htm#fixed_upper',
'sample_fixed.htm#random_opt',
'sample_fixed.htm#Positive Definite',
'sample_fixed.htm#Theory',
'sample_fixed.htm#Theory.Notation',
'sample_fixed.htm#Theory.Unconstrained Subset Covariance',
'sample_fixed.htm#Theory.Approximate Constraint Equations',
'sample_fixed.htm#Theory.Implicit Information',
'sample_fixed.htm#Theory.Implicit Covariance',
'sample_fixed.htm#Example',
'sample_fixed.htm#Other Method',
'sample_fixed.htm#2DO'
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
