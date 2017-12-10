var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'theory.htm'
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
'theory.htm#Reference',
'theory.htm#Total Likelihood',
'theory.htm#Random Likelihood, f(theta, u)',
'theory.htm#Random Likelihood, f(theta, u).Assumption',
'theory.htm#Fixed Likelihood, g(theta)',
'theory.htm#Optimal Random Effects, u^(theta)',
'theory.htm#Objective',
'theory.htm#Objective.Laplace Approximation, h(theta, u)',
'theory.htm#Objective.Laplace Objective, r(theta)',
'theory.htm#Objective.Total Objective, L(theta)',
'theory.htm#Derivative of Optimal Random Effects',
'theory.htm#Derivative of Random Constraints',
'theory.htm#Derivative of Laplace Objective',
'theory.htm#Approximate Optimal Random Effects',
'theory.htm#Approximate Optimal Random Effects.First Order, U(beta, theta, u)',
'theory.htm#Approximate Optimal Random Effects.Second Order, W(beta, theta, u)',
'theory.htm#Approximate Laplace Objective, H(beta, theta, u)',
'theory.htm#Approximate Random Constraint Function, B(beta, theta, u)',
'theory.htm#Hessian of Laplace Objective',
'theory.htm#Hessian of Random Constraints',
'theory.htm#Sparse Observed Information'
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
