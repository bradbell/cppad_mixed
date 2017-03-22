// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin initialize$$
$spell
	nr
	bool
	ldlt
	objcon
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
	rcv
	nnz
$$

$section Initialization After Constructor$$

$head Syntax$$
$icode%size_map% = %mixed_object%.initialize(%fixed_vec%, %random_vec%)%$$

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head Purpose$$
Some of the $code cppad_mixed$$ initialization requires calling the
derived class version of the
$cref ran_likelihood$$ function.
Hence this initialization cannot be done until
after the $cref/derived constructor/derived_ctor/$$ completes.

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

$head ran_likelihood_jac$$
If $cref ran_likelihood_jac$$ returns a non-empty vector,
a check is done to make sure the derivative computations are correct.

$head ran_likelihood_hes$$
If $cref ran_likelihood_hes$$ returns a non-empty vector,
a check is done to make sure the derivative computations are correct.

$head size_map$$
The return value has prototype
$codei%
	std::map<std::string, size_t> %size_map%
%$$
It represent the size of certain aspects of the problem.
For each of $code cppad_mixed$$ member variables listed below
$codei%
	%size_map%[%expression%]
%$$
where expression is one of the following:
$table
$icode expression$$ $cnext Link
$rnext
$code n_fixed_$$  $cnext
	$cref/n_fixed_/private/n_fixed_/$$
$rnext
$code n_random_$$  $cnext
	$cref/n_random_/private/n_random_/$$
$rnext
$code quasi_fixed_$$  $cnext
	$cref/quasi_fixed_/private/quasi_fixed_/$$
$rnext
$code A_rcv_.nnz()$$  $cnext
	$cref/A_rcv_/private/A_rcv_/$$
$rnext
$code A_rcv_.nr()$$  $cnext
	$cref/A_rcv_/private/A_rcv_/$$
$rnext
$code ran_like_fun_.size_var()$$  $cnext
	$cref/ran_like_fun_/private/ran_like_fun_/$$
$rnext
$code ran_like_a1fun_.size_var()$$  $cnext
	$cref/ran_like_fun_/private/ran_like_fun_/$$
$rnext
$code ran_hes_.row()$$  $cnext
	$cref/ran_hes_fun_/private/ran_hes_fun_/$$
$rnext
$code ran_hes_fun_.size_var()$$  $cnext
	$cref/ran_hes_fun_/private/ran_hes_fun_/$$
$rnext
$code 2DO$$ $cnext
	$cref/ldlt_ran_hes_/private/ldlt_ran_hes_/$$
$rnext
$code hes_cross_.row.size()$$ $cnext
	$cref/hes_cross_/private/hes_cross_/$$
$rnext
$code newton_checkpoint_.size_var()$$ $cnext
	$cref/newton_checkpoint_/private/newton_checkpoint_/$$
$rnext
$code ran_objcon_fun_.size_var()$$ $cnext
	$cref/ran_objcon_fun_/private/ran_objcon_fun_/$$
$rnext
$code ran_objcon_hes_.row.size()$$ $cnext
	$cref/ran_objcon_hes_/private/ran_objcon_fun_/ran_objcon_hes_/$$
$rnext
$code fix_like_fun.size_var()$$ $cnext
	$cref/fix_like_fun_/private/fix_like_fun_/$$
$rnext
$code fix_like_jac_.row.size()$$ $cnext
	$cref/fix_like_jac_/private/fix_like_fun_/fix_like_jac_/$$
$rnext
$code fix_like_hes_.row.size()$$ $cnext
	$cref/fix_like_hes_/private/fix_like_fun_/fix_like_hes_/$$
$rnext
$code fix_con_fun.size_var()$$ $cnext
	$cref/fix_con_fun_/private/fix_con_fun_/$$
$rnext
$code fix_con_jac_.row.size()$$ $cnext
	$cref/fix_con_jac_/private/fix_con_fun_/fix_con_jac_/$$
$rnext
$code fix_con_hes_.row.size()$$ $cnext
	$cref/fix_con_hes_/private/fix_con_fun_/fix_con_hes_/$$
$tend

$subhead num_bytes_before$$
$icode%size_map%["num_bytes_before"]%$$ is the number of bytes
allocated by $code CppAD::thread_alloc$$ and still in use
when $code initialize$$ starts.

$subhead num_bytes_after$$
$icode%size_map%["num_bytes_after"]%$$ is the number of bytes
allocated by $code CppAD::thread_alloc$$ and still in use
when $code initialization$$ is completed.

$head Example$$
The file $cref derived_ctor.cpp$$ contains an example
of using $code initialize$$.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
// ---------------------------------------------------------------------------
std::map<std::string, size_t> cppad_mixed::try_initialize(
	const d_vector&                      fixed_vec     ,
	const d_vector&                      random_vec    )
{	if( initialize_done_ )
	{	fatal_error("cppad_mixed::initialize was called twice");
	}
	if( fixed_vec.size() != n_fixed_ )
	{	fatal_error("cppad_mixed::initialize fixed_vec has wrong size");
	}
	if( random_vec.size() != n_random_ )
	{	fatal_error("cppad_mixed::initialize random_vec has wrong size");
	}
	//
	size_t thread           = CppAD::thread_alloc::thread_num();
	size_t num_bytes_before = CppAD::thread_alloc::inuse(thread);
	//
	if( n_random_ == 0 )
	{	a2_vector a2_fixed_vec( n_fixed_ ), a2_random_vec(n_random_);
		for(size_t i = 0; i < n_fixed_; i++)
			a2_fixed_vec[i] = fixed_vec[i];
		//
		a2_vector vec = ran_likelihood(a2_fixed_vec, a2_random_vec);
		if( vec.size() != 0 )
		{	std::string msg = "There are no random effects, n_random = 0,";
			msg += "\nbut ran_likelihood returns a non-empty vector";
			fatal_error(msg);
		}
	}
	else
	{
		// ran_like_
		assert( ! init_ran_like_done_ );
		init_ran_like(fixed_vec, random_vec);
		assert( init_ran_like_done_ );

		// ran_hes_
		assert( ! init_ran_hes_done_ );
		init_ran_hes(fixed_vec, random_vec);
		assert( init_ran_hes_done_ );

		// ldlt_ran_hes_
		assert( ! init_ldlt_ran_hes_done_ );
		init_ldlt_ran_hes();
		assert( init_ldlt_ran_hes_done_ );

		// hes_cross_
		assert( ! init_hes_cross_done_ );
		init_hes_cross(fixed_vec, random_vec);
		assert( init_hes_cross_done_ );

		if( ! quasi_fixed_ )
		{
			// newton_checkpoint_
			assert( ! init_newton_checkpoint_done_ );
			newton_checkpoint_.initialize(
				ran_like_a1fun_, ran_hes_rcv_, fixed_vec, random_vec
			);
			init_newton_checkpoint_done_ = true;

			// ran_objcon_fun_
			assert( ! init_ran_objcon_done_ );
			init_ran_objcon(fixed_vec, random_vec);
			assert( init_ran_objcon_done_ );

			// ran_objcon_hes_
			assert( ! init_ran_objcon_hes_done_ );
			init_ran_objcon_hes(fixed_vec, random_vec);
			assert( init_ran_objcon_hes_done_ );
		}
		// check ran_likelihood_jac
		ran_like_jac_check(fixed_vec, random_vec);
		//
		// check ran_likelihood_hes
		init_ran_hes_check(fixed_vec, random_vec);
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
	size_map["num_bytes_before"]           = num_bytes_before;
	size_map["n_fixed_"]                   = n_fixed_;
	size_map["n_random_"]                  = n_random_;
	size_map["quasi_fixed_"]               = quasi_fixed_;
	size_map["A_rcv_.nnz()"]               = A_rcv_.nnz();
	size_map["A_rcv_.nr()"]                = A_rcv_.nr();
	size_map["ran_like_fun_.size_var()"]   = ran_like_fun_.size_var();
	size_map["ran_like_a1fun_.size_var()"] = ran_like_a1fun_.size_var();
	size_map["ran_hes_rcv_.nnz()"]         = ran_hes_rcv_.nnz();
	size_map["ran_hes_fun_.size_var()"]    = ran_hes_fun_.size_var();
	size_map["hes_cross_.row.size()"]      = hes_cross_.row.size();
	size_map["newton_checkpoint_.size_var()"] = newton_checkpoint_.size_var();
	size_map["ran_objcon_fun_.size_var()"] = ran_objcon_fun_.size_var();
	size_map["ran_objcon_hes_.subset.nnz()"] = ran_objcon_hes_.subset.nnz();
	size_map["fix_like_fun_.size_var()"]   = fix_like_fun_.size_var();
	size_map["fix_like_jac_.subset.nnz()"] = fix_like_jac_.subset.nnz();
	size_map["fix_like_hes_.subset.nnz()"] = fix_like_hes_.subset.nnz();
	size_map["fix_con_fun_.size_var()"]    = fix_con_fun_.size_var();
	size_map["fix_con_jac_.subset.nnz()"]    = fix_con_jac_.subset.nnz();
	size_map["fix_con_hes_.subset.nnz()"]    = fix_con_hes_.subset.nnz();
	size_map["num_bytes_after"]            = CppAD::thread_alloc::inuse(thread);
	return size_map;
}
// ---------------------------------------------------------------------------
std::map<std::string, size_t> cppad_mixed::initialize(
	const d_vector&                       fixed_vec      ,
	const d_vector&                       random_vec     )
{	std::map<std::string, size_t> ret;
	try
	{	ret = try_initialize(fixed_vec, random_vec);
	}
	catch(const CppAD::mixed::exception& e)
	{	std::string error_message = e.message("initialize");
		fatal_error(error_message);
		assert(false);
	}
	return ret;
}
