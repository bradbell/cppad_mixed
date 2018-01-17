var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'install_unix.htm',
'run_cmake.sh.htm'
];
var list_down2 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_18.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down1 = [
'example_install.sh.htm',
'run_cmake.sh.htm',
'check_install.sh.htm'
];
var list_current0 = [
'run_cmake.sh.htm#verbose_makefile',
'run_cmake.sh.htm#build_type',
'run_cmake.sh.htm#Prefixes',
'run_cmake.sh.htm#Prefixes.Debug and Release',
'run_cmake.sh.htm#cppad_cxx_flags',
'run_cmake.sh.htm#cmake_libdir',
'run_cmake.sh.htm#ldlt_cholmod',
'run_cmake.sh.htm#use_atomic_cholesky',
'run_cmake.sh.htm#checkpoint_newton_step',
'run_cmake.sh.htm#optimize_cppad_function',
'run_cmake.sh.htm#hide_ipopt_scaling',
'run_cmake.sh.htm#for_hes_sparsity',
'run_cmake.sh.htm#Testing Speed and Memory'
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
