// Child table for section ipopt_xam
document.write('\
<select onchange="ipopt_xam_child(this)">\
<option>ipopt_xam-&gt;</option>\
<option>ipopt_nlp_xam</option>\
<option>ipopt_xam_ctor</option>\
<option>ipopt_xam_get_nlp_info</option>\
<option>ipopt_xam_get_bounds_info</option>\
<option>ipopt_xam_get_starting_point</option>\
<option>ipopt_xam_eval_f</option>\
<option>ipopt_xam_eval_grad_f</option>\
<option>ipopt_xam_eval_g</option>\
<option>ipopt_xam_eval_jac_g</option>\
<option>ipopt_xam_eval_h</option>\
<option>ipopt_xam_finalize_solution</option>\
<option>ipopt_xam_intermediate_callback</option>\
<option>ipopt_run_xam</option>\
</select>\
');
function ipopt_xam_child(item)
{	var child_list = [
		'ipopt_nlp_xam.htm',
		'ipopt_xam_ctor.htm',
		'ipopt_xam_get_nlp_info.htm',
		'ipopt_xam_get_bounds_info.htm',
		'ipopt_xam_get_starting_point.htm',
		'ipopt_xam_eval_f.htm',
		'ipopt_xam_eval_grad_f.htm',
		'ipopt_xam_eval_g.htm',
		'ipopt_xam_eval_jac_g.htm',
		'ipopt_xam_eval_h.htm',
		'ipopt_xam_finalize_solution.htm',
		'ipopt_xam_intermediate_callback.htm',
		'ipopt_run_xam.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
