var list_across0 = [
'_contents.htm',
'_reference.htm',
'_index.htm',
'_search.htm',
'_external.htm'
];
var list_up0 = [
'cppad_mixed.htm',
'whats_new_17.htm'
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
var list_down0 = [
'whats_new_15.htm',
'whats_new_16.htm'
];
var list_current0 = [
'whats_new_17.htm#Contents',
'whats_new_17.htm#12-10',
'whats_new_17.htm#10-27',
'whats_new_17.htm#10-24',
'whats_new_17.htm#10-09',
'whats_new_17.htm#10-07',
'whats_new_17.htm#09-30',
'whats_new_17.htm#09-23',
'whats_new_17.htm#09-21',
'whats_new_17.htm#09-18',
'whats_new_17.htm#09-16',
'whats_new_17.htm#09-15',
'whats_new_17.htm#09-14',
'whats_new_17.htm#09-02',
'whats_new_17.htm#08-30',
'whats_new_17.htm#08-01',
'whats_new_17.htm#04-24',
'whats_new_17.htm#04-23',
'whats_new_17.htm#04-06',
'whats_new_17.htm#04-02',
'whats_new_17.htm#03-27',
'whats_new_17.htm#03-25',
'whats_new_17.htm#03-23',
'whats_new_17.htm#03-20',
'whats_new_17.htm#03-12',
'whats_new_17.htm#03-11',
'whats_new_17.htm#03-10',
'whats_new_17.htm#03-09',
'whats_new_17.htm#03-08',
'whats_new_17.htm#03-06',
'whats_new_17.htm#03-02',
'whats_new_17.htm#03-01',
'whats_new_17.htm#01-26',
'whats_new_17.htm#01-24',
'whats_new_17.htm#01-22',
'whats_new_17.htm#01-14'
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
