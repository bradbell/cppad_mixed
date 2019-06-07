// Child table for section ipopt_fixed
document.write('\
<select onchange="ipopt_fixed_child(this)">\
<option>ipopt_fixed-&gt;</option>\
<option>ipopt_fixed_ctor</option>\
<option>ipopt_fixed_get_nlp_info</option>\
<option>ipopt_fixed_get_bounds_info</option>\
<option>ipopt_fixed_get_starting_point</option>\
<option>ipopt_fixed_eval_f</option>\
<option>ipopt_fixed_eval_grad_f</option>\
<option>ipopt_fixed_eval_g</option>\
<option>ipopt_fixed_eval_jac_g</option>\
<option>ipopt_fixed_eval_h</option>\
<option>ipopt_fixed_finalize_solution</option>\
<option>ipopt_fixed_adaptive_derivative_check</option>\
<option>ipopt_fixed_new_random</option>\
<option>ipopt_xam</option>\
</select>\
');
function ipopt_fixed_child(item)
{	var child_list = [
		'ipopt_fixed_ctor.htm',
		'ipopt_fixed_get_nlp_info.htm',
		'ipopt_fixed_get_bounds_info.htm',
		'ipopt_fixed_get_starting_point.htm',
		'ipopt_fixed_eval_f.htm',
		'ipopt_fixed_eval_grad_f.htm',
		'ipopt_fixed_eval_g.htm',
		'ipopt_fixed_eval_jac_g.htm',
		'ipopt_fixed_eval_h.htm',
		'ipopt_fixed_finalize_solution.htm',
		'ipopt_fixed_adaptive_derivative_check.htm',
		'ipopt_fixed_new_random.htm',
		'ipopt_xam.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
