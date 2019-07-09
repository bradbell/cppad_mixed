// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-19 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin initialize$$
$spell
	iterator
	itr
	cout
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
This $code cppad_mixed$$ $cref base_class$$ member function is public.

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

$head size_map$$
The return value has prototype
$codei%
	std::map<std::string, size_t> %size_map%
%$$
It represent the size of certain aspects of the problem.
The actual fields in $icode size_map$$ are not specified,
but they can be inspected. For example,
$codei%
	std::map<std::string, size_t>::iterator itr;
	for(itr = %size_map%.begin(); itr != %size_map%.end(); itr++)
		std::cout << itr->first << " = " itr->second << "\n";
%$$

$head init_newton_checkpoint_done_$$
The input value of this member variable must be false.
If $code quasi_fixed_$$ is true, its output value is false,
otherwise its output value is true.

$head Example$$
The file $cref derived_ctor.cpp$$ contains an example
of using $code initialize$$.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
// ---------------------------------------------------------------------------
std::map<std::string, size_t> cppad_mixed::try_initialize(
	const d_vector&                      fixed_ini     ,
	const d_vector&                      random_ini    )
{	if( initialize_done_ )
	{	fatal_error("cppad_mixed::initialize was called twice");
	}
	if( fixed_ini.size() != n_fixed_ )
	{	fatal_error("cppad_mixed::initialize fixed_vec has wrong size");
	}
	if( random_ini.size() != n_random_ )
	{	fatal_error("cppad_mixed::initialize random_vec has wrong size");
	}
	// convert fixed_ini -> fixed_vec and random_ini -> random_vec
	//
	// Eigen will short curcit operations when values are zero. This can
	// result in a different operation sequence. Create versions of the
	// input vectors with no zeros in them.
	double tiny = std::numeric_limits<double>::min();
	assert( tiny > 0.0 );
	d_vector fixed_vec  = fixed_ini;
	d_vector random_vec = random_ini;
	for(size_t i = 0; i < n_fixed_; i++)
	{	if( fixed_vec[i] == 0.0 )
			fixed_vec[i] = tiny;
	}
	for(size_t i = 0; i < n_random_; i++)
	{	if( random_vec[i] == 0.0 )
			random_vec[i] = tiny;
	}
	// ------------------------------------------------------------------------
	//
	size_t thread           = CppAD::thread_alloc::thread_num();
	size_t num_bytes_before = CppAD::thread_alloc::inuse(thread);
	//
	if( n_random_ == 0 )
	{	a1_vector a1_fixed_vec( n_fixed_ ), a1_random_vec(n_random_);
		for(size_t i = 0; i < n_fixed_; i++)
			a1_fixed_vec[i] = fixed_vec[i];
		for(size_t i = 0; i < n_random_; i++)
			a1_random_vec[i] = random_vec[i];
		//
		a1_vector vec = ran_likelihood(a1_fixed_vec, a1_random_vec);
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

		// ran_jac_
		assert( ! init_ran_jac_done_ );
		init_ran_jac(fixed_vec, random_vec);
		assert( init_ran_jac_done_ );

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

		assert( ! init_laplace_obj_done_ );
		if( ! quasi_fixed_ )
		{
			// laplace_obj_fun_
			assert( ! init_laplace_obj_done_ );
			init_laplace_obj(fixed_vec, random_vec);
			assert( init_laplace_obj_done_ );

			// laplace_obj_hes_
			assert( ! init_laplace_obj_hes_done_ );
			init_laplace_obj_hes(fixed_vec, random_vec);
			assert( init_laplace_obj_hes_done_ );
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
	size_map["num_bytes_before"]           = num_bytes_before;
	size_map["n_fixed_"]                   = n_fixed_;
	size_map["n_random_"]                  = n_random_;
	size_map["quasi_fixed_"]               = quasi_fixed_;
	size_map["A_rcv_.nnz()"]               = A_rcv_.nnz();
	size_map["A_rcv_.nr()"]                = A_rcv_.nr();
	size_map["ran_like_fun_.size_var()"]   = ran_like_fun_.size_var();
	size_map["ran_like_a1fun_.size_var()"] = ran_like_a1fun_.size_var();
	size_map["ran_hes_uu_rcv_.nnz()"]      = ran_hes_uu_rcv_.nnz();
	size_map["ran_hes_fun_.size_var()"]    = ran_hes_fun_.size_var();
	size_map["hes_cross_.subset.nnz()"]    = hes_cross_.subset.nnz();
	size_map["laplace_obj_fun_.size_var()"] = laplace_obj_fun_.size_var();
	size_map["laplace_obj_hes_.subset.nnz()"] = laplace_obj_hes_.subset.nnz();
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
# ifndef NDEBUG
	ret = try_initialize(fixed_vec, random_vec);
# else
	try
	{	ret = try_initialize(fixed_vec, random_vec);
	}
	catch(const CppAD::mixed::exception& e)
	{	std::string error_message = e.message("initialize");
		fatal_error(error_message);
		assert(false);
	}
# endif
	return ret;
}
