// Child table for section ldlt_cholmod
document.write('\
<select onchange="ldlt_cholmod_child(this)">\
<option>ldlt_cholmod-&gt;</option>\
<option>ldlt_cholmod_ctor</option>\
<option>ldlt_cholmod_dtor</option>\
<option>ldlt_cholmod_init</option>\
<option>ldlt_cholmod_pattern</option>\
<option>ldlt_cholmod_update</option>\
<option>ldlt_cholmod_logdet</option>\
<option>ldlt_cholmod_solve_H</option>\
<option>ldlt_cholmod_sim_cov</option>\
<option>ldlt_cholmod_inv</option>\
<option>ldlt_cholmod.cpp</option>\
<option>cholmod_solve_xam</option>\
<option>cholmod_solve2_a.cpp</option>\
<option>cholmod_solve2_sim.cpp</option>\
</select>\
');
function ldlt_cholmod_child(item)
{	var child_list = [
		'ldlt_cholmod_ctor.htm',
		'ldlt_cholmod_dtor.htm',
		'ldlt_cholmod_init.htm',
		'ldlt_cholmod_pattern.htm',
		'ldlt_cholmod_update.htm',
		'ldlt_cholmod_logdet.htm',
		'ldlt_cholmod_solve_h.htm',
		'ldlt_cholmod_sim_cov.htm',
		'ldlt_cholmod_inv.htm',
		'ldlt_cholmod.cpp.htm',
		'cholmod_solve_xam.htm',
		'cholmod_solve2_a.cpp.htm',
		'cholmod_solve2_sim.cpp.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
