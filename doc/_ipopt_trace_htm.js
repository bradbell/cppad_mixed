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
'optimize_fixed.htm',
'ipopt_trace.htm'
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
'public.htm',
'private.htm'
];
var list_down2 = [
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
var list_down1 = [
'optimize_fixed.cpp.htm',
'ipopt_options.htm',
'ipopt_trace.htm'
];
var list_current0 = [
'ipopt_trace.htm#iter',
'ipopt_trace.htm#objective',
'ipopt_trace.htm#inf_pr',
'ipopt_trace.htm#inf_du',
'ipopt_trace.htm#lg(mu)',
'ipopt_trace.htm#||d||',
'ipopt_trace.htm#lg(rg)',
'ipopt_trace.htm#alpha_du',
'ipopt_trace.htm#alpha_pr',
'ipopt_trace.htm#ls',
'ipopt_trace.htm#Reference'
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