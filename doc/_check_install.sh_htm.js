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
'check_install.sh.htm'
];
var list_down2 = [
'install_unix.htm',
'theory.htm',
'base_class.htm',
'namespace.htm',
'user.htm',
'whats_new_17.htm',
'wish_list.htm',
'math_notation.htm'
];
var list_down1 = [
'example_install.sh.htm',
'run_cmake.sh.htm',
'check_install.sh.htm'
];
var list_current0 = [
'check_install.sh.htm#Syntax',
'check_install.sh.htm#build_type',
'check_install.sh.htm#Prefixes',
'check_install.sh.htm#cmake_libdir',
'check_install.sh.htm#example_file',
'check_install.sh.htm#PKG_CONFIG_PATH',
'check_install.sh.htm#LD_LIBRARY_PATH',
'check_install.sh.htm#Create Temporary',
'check_install.sh.htm#example_name',
'check_install.sh.htm#main',
'check_install.sh.htm#gsl_libs',
'check_install.sh.htm#ipopt_libs',
'check_install.sh.htm#suitesparse_libs',
'check_install.sh.htm#Compile and Link',
'check_install.sh.htm#Run Example'
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
