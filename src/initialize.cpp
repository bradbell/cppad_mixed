// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>

/*
$begin initialize$$
$spell
	CppAD
	ranobj
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

$subhead hes_ran_row$$
$icode%size_map%["hes_ran_row"]%$$ is the
number of non-zero entries in hessian of
random likelihood w.r.t random effects
$latex f_{uu}^{(2)} ( \theta , u )$$.

$subhead hes_cross_row$$
$icode%size_map%["hes_cross_row"]%$$ is the
number of non-zero cross entries in hessian of
random likelihood w.r.t fixed and random effects
$latex f_{u \theta}^{(2)} ( \theta , u )$$.

$subhead ran_likelihood_fun$$
$icode%size_map%["ran_likelihood_fun"]%$$ is the
number of variables in the $code ADFun<double>$$ version of
$cref ran_likelihood$$.

$subhead ran_likelihood_a1fun$$
$icode%size_map%["ran_likelihood_a1fun"]%$$ is the
number of variables in the $code ADFun<a1_double>$$ version of
$cref ran_likelihood$$.

$subhead hes_ran_fun$$
$icode%size_map%["hes_ran_fun_"]%$$ is the
number of variables in the $code ADFun<double>$$ object that
computes the Hessian with respect to the random effects.

$subhead newton_atom$$
$icode%size_map%["newton_atom"]%$$ is the
number of variables in the $code ADFun<double>$$
object used to evaluate the newton step
$latex \[
	s = f_{uu}^{(2)} ( \theta , u )^{-1} v
\] $$
and the log of the determinant of the matrix being inverted.

$subhead ranobj_fun$$
$icode%size_map%["ranobj_2"]%$$ is the
number of variables in the $code ADFun<double>$$
object used to evaluate Hessian, w.r.t. fixed effects, of the objective
$latex r^{(2)} ( \theta )$$.

$subhead fix_likelihood_fun$$
$icode%size_map%["fix_likelihood"]%$$ is the
number of variables in the $code ADFun<double>$$ function
that computes the fixed likelihood $latex g( \theta )$$.

$subhead constraint$$
$icode%size_map%["constraint"]%$$ is the
number of variables in the $code ADFun<double>$$ function
that computes the general constraints for the fixed effects.

$subhead num_bytes_before$$
$icode%size_map%["num_bytes_before"]%$$ is the
number of bytes allocated by
$code CppAD::thread_alloc$$ before the initialization is started.

$subhead num_bytes_after$$
$icode%size_map%["num_bytes_after"]%$$ is the
number of bytes allocated by
$code CppAD::thread_alloc$$ after the initialization is completed.

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
		assert( ! init_ran_likelihood_done_ );
		init_ran_likelihood(fixed_vec, random_vec);
		assert( init_ran_likelihood_done_ );

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
			newton_atom_.initialize(ran_likelihood_a1fun_, fixed_vec, random_vec);
			record_newton_atom_done_ = true;

			// ranobj_fun_
			assert( ! init_ranobj_done_ );
			init_ranobj(fixed_vec, random_vec);
			assert( init_ranobj_done_ );

			// hes_ranobj_
			assert( ! init_hes_ranobj_done_ );
			init_hes_ranobj(fixed_vec, random_vec);
			assert( init_hes_ranobj_done_ );
		}
	}

	// fix_likelihood_fun_
	assert( ! init_fix_likelihood_done_ );
	init_fix_likelihood(fixed_vec);
	assert( init_fix_likelihood_done_ );

	// fix_constraint_fun_
	assert( ! init_fix_constraint_done_ );
	init_fix_constraint(fixed_vec);
	assert( init_fix_constraint_done_ );

	// initialize_done_
	initialize_done_ = true;

	// return value
	std::map<std::string, size_t> size_map;
	size_map["hes_ran_row"]    = hes_ran_.row.size();
	size_map["hes_cross_row"]  = hes_cross_.row.size();
	size_map["ran_likelihood_fun"]   = ran_likelihood_fun_.size_var();
	size_map["ran_likelihood_a1fun"] = ran_likelihood_a1fun_.size_var();
	size_map["hes_ran_fun"]    = hes_ran_fun_.size_var();
	size_map["newton_atom"]    = newton_atom_.size_var();
	size_map["ranobj_fun"]    = ranobj_fun_.size_var();
	size_map["fix_likelihood_fun"]   = fix_likelihood_fun_.size_var();
	size_map["constraint"]     = fix_constraint_fun_.size_var();
	//
	size_map["num_bytes_before"]  = num_bytes_before;
	size_map["num_bytes_after"]   = CppAD::thread_alloc::inuse(thread);
	return size_map;
}

