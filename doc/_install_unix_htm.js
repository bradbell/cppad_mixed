var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'install_unix.htm'
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
var list_down0 = [
'example_install.sh.htm',
'run_cmake.sh.htm',
'check_install.sh.htm'
];
var list_current0 = [
'install_unix.htm#System Requirements',
'install_unix.htm#System Requirements.C++ Compiler',
'install_unix.htm#System Requirements.git',
'install_unix.htm#System Requirements.cmake',
'install_unix.htm#System Requirements.pkg-config',
'install_unix.htm#System Requirements.wget',
'install_unix.htm#System Requirements.Fortran Compiler',
'install_unix.htm#System Requirements.gsl',
'install_unix.htm#System Requirements.suitesparse',
'install_unix.htm#Download',
'install_unix.htm#Paths',
'install_unix.htm#Paths.PKG_CONFIG_PATH',
'install_unix.htm#Paths.LD_LIBRARY_PATH',
'install_unix.htm#Special Requirements',
'install_unix.htm#Special Requirements.run_cmake.sh',
'install_unix.htm#Special Requirements.eigen',
'install_unix.htm#Special Requirements.Ipopt',
'install_unix.htm#Special Requirements.CppAD',
'install_unix.htm#cppad_mixed',
'install_unix.htm#cppad_mixed.Cmake Command',
'install_unix.htm#cppad_mixed.Check',
'install_unix.htm#cppad_mixed.Speed',
'install_unix.htm#cppad_mixed.Install',
'install_unix.htm#Example',
'install_unix.htm#Example.Installation',
'install_unix.htm#Example.Linking',
'install_unix.htm#Example.Using cppad_mixed'
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
