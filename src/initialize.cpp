// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>

/*
$begin initialize$$
$spell
	init
	CppAD
	ran_obj
	logdet
	cppad
	std
	obj
	vec
	const
	Cpp
	var
	hes
	num
	alloc
	jac
	Jacobian
$$

$section Initialization After Constructor$$

$head Syntax$$
$icode%size_map% = %mixed_object%.initialize(%fixed_vec%, %random_vec%)%$$

$head Purpose$$
Some of the $code cppad_mixed$$ initialization requires calling the
derived class version of the
$cref ran_likelihood$$ function.
Hence this initialization cannot be done until
after the $cref/derived constructor/derived_ctor/$$ completes.

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which certain $code CppAD::ADFun$$
objects are recorded.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which certain $code CppAD::ADFun$$
objects are recorded.

$head size_map$$
The return value has prototype
$codei%
	std::map<std::string, size_t> %size_map%
%$$
It represent the size of certain aspects of the problem.
For each of $code cppad_mixed$$ member variables listed below
$codei%
	%size_map%[%member_variable%]
%$$
is the corresponding size.

$subhead fix_con_fun_$$
$icode%size_map%["fix_con_fun_"]%$$ is the number of variables in
$cref/fix_con_fun_/init_fix_con/fix_con_fun_/$$.

$subhead fix_con_hes_$$
$icode%size_map%["fix_con_hes_"]%$$ is the size of the row vector in
$cref/fix_con_hes_/init_fix_con/fix_con_hes_/$$.
This is also the number of non-zeros
in the Hessian $latex c_{\theta,\theta} ( \theta )$$.

$subhead fix_con_jac_$$
$icode%size_map%["fix_con_jac_"]%$$ is the size of the row vector in
$cref/fix_con_jac_/init_fix_con/fix_con_jac_/$$.
This is also the number of non-zeros
in the Jacobian $latex c_\theta ( \theta )$$.

$subhead fix_like_fun_$$
$icode%size_map%["fix_like_fun_"]%$$ is the number of variables in
$cref/fix_like_fun_/init_fix_like/fix_like_fun_/$$.

$subhead fix_like_hes_$$
$icode%size_map%["fix_like_hes_"]%$$ is the size of the row vector in
$cref/fix_like_hes_/init_fix_like/fix_like_hes_/$$.
This is also the number of non-zeros
in the Hessian $latex g_{\theta,\theta} ( \theta )$$.

$subhead fix_like_jac_$$
$icode%size_map%["fix_like_jac_"]%$$ is the size of the row vector in
$cref/fix_like_jac_/init_fix_like/fix_like_jac_/$$.
This is also the number of non-zeros
in the Jacobian $latex g_\theta ( \theta )$$.

$subhead hes_ran_$$
$icode%size_map%["hes_ran_"]%$$ is the size of the row vector in
$cref/hes_ran_/init_hes_ran/hes_ran_/$$.
This is also the number of non-zeros
in the lower triangle of the Hessian
$latex \[
	f_{u,u} ( \theta , u )
\]$$
see $cref/f(theta, u)/
	theory/
	Random Likelihood, f(theta, u)
/$$

$subhead hes_cross_$$
$icode%size_map%["hes_cross_"]%$$ is the size of the row vector in
$cref/hes_cross_/init_hes_cross/hes_cross_/$$.
This is also the number of non-zeros
in the Hessian
$latex \[
	f_{u,\theta} ( \theta , u )
\]$$
see $cref/f(theta, u)/
	theory/
	Random Likelihood, f(theta, u)
/$$

$subhead hes_ran_fun_$$
$icode%size_map%["hes_ran_fun_"]%$$ is the number of variables in
$cref/hes_ran_fun_/init_hes_ran/hes_ran_fun_/$$.

$subhead hes_ran_obj_$$
$icode%size_map%["hes_ran_obj_"]%$$ is the size of the row vector in
$cref/hes_ran_obj_/init_hes_ran_obj/hes_ran_obj_/$$.
This is also the number of non-zeros
in the lower triangle of the Hessian
$latex \[
	r_{\theta,\theta} ( \theta )
	=
	H_{\beta,\beta} [ \beta, \theta , \hat{u} ( \theta) ]
\] $$
see
$cref/H(beta, theta, u)
	/theory
	/Approximate Random Objective, H(beta, theta, u)
/$$.

$subhead newton_atom_$$
If $cref/quasi_fixed_/private/quasi_fixed_/$$ is false,
$icode%size_map%["newton_atom_"]%$$ is the number of variables in the
$cref/newton_atom_/private/newton_atom_/$$
tape used to represent the $cref newton_step$$ checkpoint function.

$subhead num_bytes_after$$
$icode%size_map%["num_bytes_after"]%$$ is the number of bytes
allocated by $code CppAD::thread_alloc$$ and still in use
when $code initialization$$ is completed.

$subhead num_bytes_before$$
$icode%size_map%["num_bytes_before"]%$$ is the number of bytes
allocated by $code CppAD::thread_alloc$$ and still in use
when $code initialize$$ starts.

$subhead ran_like_a1fun_$$
$icode%size_map%["ran_like_a1fun_"]%$$ is the number of variables in
$cref/ran_like_a1fun_/init_ran_like/ran_like_a1fun_/$$.

$subhead ran_like_fun_$$
$icode%size_map%["ran_like_fun_"]%$$ is the number of variables in
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$.

$subhead ran_obj_fun_$$
$icode%size_map%["ran_obj_fun_"]%$$ is the number of variables in
$cref/ran_obj_fun_/init_ran_obj/ran_obj_fun_/$$.


$head Example$$
The file $cref derived_xam.cpp$$ contains an example
of using $code initialize$$.

$end
*/


std::map<std::string, size_t> cppad_mixed::initialize(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	if( initialize_done_ )
	{	fatal_error("cppad_mixed::initialize was called twice");
	}
	size_t thread           = CppAD::thread_alloc::thread_num();
	size_t num_bytes_before = CppAD::thread_alloc::inuse(thread);
	//
	if( n_random_ == 0 )
	{	a1d_vector a1_fixed_vec( n_fixed_ ), a1_random_vec;
		for(size_t i = 0; i < n_fixed_; i++)
			a1_fixed_vec[i] = fixed_vec[i];
		a1d_vector vec = ran_likelihood(a1_fixed_vec, a1_random_vec);
		if( vec.size() != 0 )
		{	std::string msg = "There are no random effects (n_random = 0),";
			msg += "\nbut the function ran_likelihood returns a non-empty vector";
			fatal_error(msg);
		}
	}
	else
	{
		// ran_likelihood_
		assert( ! init_ran_like_done_ );
		init_ran_like(fixed_vec, random_vec);
		assert( init_ran_like_done_ );

		// hes_ran_
		assert( ! init_hes_ran_done_ );
		init_hes_ran(fixed_vec, random_vec);
		assert( init_hes_ran_done_ );

		// hes_cross_
		assert( ! init_hes_cross_done_ );
		init_hes_cross(fixed_vec, random_vec);
		assert( init_hes_cross_done_ );

		if( ! quasi_fixed_ )
		{
			// newton_atom_
			assert( ! record_newton_atom_done_ );
			newton_atom_.initialize(
				ran_like_a1fun_, fixed_vec, random_vec
			);
			record_newton_atom_done_ = true;

			// ran_obj_fun_
			assert( ! init_ran_obj_done_ );
			init_ran_obj(fixed_vec, random_vec);
			assert( init_ran_obj_done_ );

			// hes_ran_obj_
			assert( ! init_hes_ran_obj_done_ );
			init_hes_ran_obj(fixed_vec, random_vec);
			assert( init_hes_ran_obj_done_ );
		}
	}

	// fix_like_fun_
	assert( ! init_fix_like_done_ );
	init_fix_like(fixed_vec);
	assert( init_fix_like_done_ );

	// fix_con_fun_
	assert( ! init_fix_con_done_ );
	init_fix_con(fixed_vec);
	assert( init_fix_con_done_ );

	// initialize_done_
	initialize_done_ = true;

	// return value
	std::map<std::string, size_t> size_map;
	size_map["num_bytes_before"]      = num_bytes_before;
	size_map["ran_like_fun_"]   = ran_like_fun_.size_var();
	size_map["ran_like_a1fun_"] = ran_like_a1fun_.size_var();
	//
	size_map["hes_ran_"]              = hes_ran_.row.size();
	size_map["hes_ran_fun_"]          = hes_ran_fun_.size_var();
	//
	size_map["hes_cross_"]            = hes_cross_.row.size();
	if( ! quasi_fixed_ )
	{	size_map["newton_step_"]      = newton_atom_.size_var();
		//
		size_map["ran_obj_fun_"]       = ran_obj_fun_.size_var();
		size_map["hes_ran_obj_"]       = hes_ran_obj_.row.size();
	}
	size_map["fix_like_fun_"]   = fix_like_fun_.size_var();
	size_map["fix_like_jac_"]   = fix_like_jac_.row.size();
	size_map["fix_like_hes_"]   = fix_like_hes_.row.size();
	//
	size_map["fix_con_fun_"]   = fix_con_fun_.size_var();
	size_map["fix_con_jac_"]   = fix_con_jac_.row.size();
	size_map["fix_con_hes_"]   = fix_con_hes_.row.size();
	//
	size_map["num_bytes_after"]       = CppAD::thread_alloc::inuse(thread);
	return size_map;
}

