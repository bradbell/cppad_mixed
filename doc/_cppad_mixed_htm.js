var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm'
];
var list_down0 = [
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
'cppad_mixed.htm#License',
'cppad_mixed.htm#Source Code Repository',
'cppad_mixed.htm#Notation',
'cppad_mixed.htm#Notation.Fixed Effects, theta',
'cppad_mixed.htm#Notation.Random Effects, u',
'cppad_mixed.htm#Notation.Data, y, z',
'cppad_mixed.htm#Notation.Fixed Prior Density, p(theta)',
'cppad_mixed.htm#Notation.Fixed Data Density, p(z|theta)',
'cppad_mixed.htm#Notation.Random Prior Density, p(u|theta)',
'cppad_mixed.htm#Notation.Random Data Density, p(y|theta,u)',
'cppad_mixed.htm#Notation.Fixed Constraint Function, c(theta)',
'cppad_mixed.htm#Notation.Optimal Random Effects, u^(theta)',
'cppad_mixed.htm#Notation.Random Constraint Matrix, A',
'cppad_mixed.htm#Notation.Random Constraint Function, A*u^(theta)',
'cppad_mixed.htm#Problem',
'cppad_mixed.htm#Problem.Maximum Likelihood',
'cppad_mixed.htm#Problem.No Random Effects',
'cppad_mixed.htm#Problem.Fixed Constraints, c',
'cppad_mixed.htm#Problem.Random Constraints',
'cppad_mixed.htm#Negative Log-Density Vector',
'cppad_mixed.htm#Contents'
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
