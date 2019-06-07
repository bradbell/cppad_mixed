// Child table for section ipopt_random
document.write('\
<select onchange="ipopt_random_child(this)">\
<option>ipopt_random-&gt;</option>\
<option>ipopt_random_ctor</option>\
<option>ipopt_random_get_nlp_info</option>\
<option>ipopt_random_get_bounds_info</option>\
<option>ipopt_random_get_starting_point</option>\
<option>ipopt_random_eval_f</option>\
<option>ipopt_random_eval_grad_f</option>\
<option>ipopt_random_eval_g</option>\
<option>ipopt_random_eval_jac_g</option>\
<option>ipopt_random_eval_h</option>\
<option>ipopt_random_finalize_solution</option>\
</select>\
');
function ipopt_random_child(item)
{	var child_list = [
		'ipopt_random_ctor.htm',
		'ipopt_random_get_nlp_info.htm',
		'ipopt_random_get_bounds_info.htm',
		'ipopt_random_get_starting_point.htm',
		'ipopt_random_eval_f.htm',
		'ipopt_random_eval_grad_f.htm',
		'ipopt_random_eval_g.htm',
		'ipopt_random_eval_jac_g.htm',
		'ipopt_random_eval_h.htm',
		'ipopt_random_finalize_solution.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
