// Child table for section private
document.write('\
<select onchange="private_child(this)">\
<option>private-&gt;</option>\
<option>pack</option>\
<option>unpack</option>\
<option>init_ran_jac</option>\
<option>init_ran_hes</option>\
<option>init_laplace_obj</option>\
<option>init_ldlt_ran_hes</option>\
<option>init_fix_con</option>\
<option>init_fix_like</option>\
<option>init_hes_cross</option>\
<option>init_laplace_obj_hes</option>\
<option>init_ran_like</option>\
<option>fix_con_eval</option>\
<option>fix_con_hes</option>\
<option>fix_con_jac</option>\
<option>fix_like_eval</option>\
<option>fix_like_hes</option>\
<option>fix_like_jac</option>\
<option>logdet_jac</option>\
<option>ran_like_hes</option>\
<option>ran_con_eval</option>\
<option>ran_con_jac</option>\
<option>ran_obj_eval</option>\
<option>ran_obj_jac</option>\
<option>laplace_obj_hes</option>\
<option>update_factor</option>\
</select>\
');
function private_child(item)
{	var child_list = [
		'pack.htm',
		'unpack.htm',
		'init_ran_jac.htm',
		'init_ran_hes.htm',
		'init_laplace_obj.htm',
		'init_ldlt_ran_hes.htm',
		'init_fix_con.htm',
		'init_fix_like.htm',
		'init_hes_cross.htm',
		'init_laplace_obj_hes.htm',
		'init_ran_like.htm',
		'fix_con_eval.htm',
		'fix_con_hes.htm',
		'fix_con_jac.htm',
		'fix_like_eval.htm',
		'fix_like_hes.htm',
		'fix_like_jac.htm',
		'logdet_jac.htm',
		'ran_like_hes.htm',
		'ran_con_eval.htm',
		'ran_con_jac.htm',
		'ran_obj_eval.htm',
		'ran_obj_jac.htm',
		'laplace_obj_hes.htm',
		'update_factor.htm'
	];
	var index = item.selectedIndex;
	item.selectedIndex = 0;
	if(index > 0)
		document.location = child_list[index-1];
}
