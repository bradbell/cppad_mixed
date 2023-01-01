// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
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

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
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

$subhead n_fixed$$
the number of fixed effects.

$subhead n_random$$
the number of fixed effects.

$subhead quasi_fixed$$
If this is one (zero) are a using a quasi-Newton (Newton) method
for optimizing the fixed effects.

$subhead A_nr$$
is the number of rows in the liner constraint matrix A
(the matrix has $icode n_fixed$$ columns).

$subhead A_nnz$$
is the number of non-zeros in the liner constraint matrix A.

$subhead ran_like_fun.size_var$$
is the number of variables in the algorithm that maps the
fixed and random effects to the part of the likelihood that depend
on the random effects.

$subhead fix_like_fun.size_var$$
is the number of variables in the algorithm that maps the
fixed effects to the part of the likelihood that does not depend
on the random effects.

$subhead Other Fields$$
Not all the fields in $icode size_map$$ are specified,
but they can be inspected. For example,
$codei%
	std::map<std::string, size_t>::iterator itr;
	for(itr = %size_map%.begin(); itr != %size_map%.end(); itr++)
		std::cout << itr->first << " = " << itr->second << "\n";
%$$

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
{
	if( trace_init_ )
		std::cout << "Begin cppad_mixed::initialize\n";
	//
	if( initialize_done_ )
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
		if( trace_init_ )
			std::cout << "No random effects\n";
	}
	else
	{
		// ran_like_
		assert( ! init_ran_like_done_ );
		init_ran_like(fixed_vec, random_vec);
		assert( init_ran_like_done_ );
		if( trace_init_ )
			std::cout << "init_ran_like_done_\n";

		// ran_jac_
		assert( ! init_ran_jac_done_ );
		init_ran_jac(fixed_vec, random_vec);
		assert( init_ran_jac_done_ );
		if( trace_init_ )
			std::cout << "init_ran_jac_done_\n";

		// ran_hes_
		assert( ! init_ran_hes_done_ );
		init_ran_hes(fixed_vec, random_vec);
		assert( init_ran_hes_done_ );
		if( trace_init_ )
			std::cout << "init_ran_hes_done_\n";

		// ldlt_ran_hes_
		assert( ! init_ldlt_ran_hes_done_ );
		init_ldlt_ran_hes();
		assert( init_ldlt_ran_hes_done_ );
		if( trace_init_ )
			std::cout << "init_ldlt_ran_hes_done_\n";

		// hes_cross_
		assert( ! init_hes_cross_done_ );
		init_hes_cross(fixed_vec, random_vec);
		assert( init_hes_cross_done_ );
		if( trace_init_ )
			std::cout << "init_hes_cross_done_\n";
	}

	// fix_like_fun_
	assert( ! init_fix_like_done_ );
	init_fix_like(fixed_vec);
	assert( init_fix_like_done_ );
	if( trace_init_ )
		std::cout << "init_fix_like_done_\n";

	// fix_con_fun_
	assert( ! init_fix_con_done_ );
	init_fix_con(fixed_vec);
	assert( init_fix_con_done_ );
	if( trace_init_ )
		std::cout << "init_fix_con_done_\n";

	// initialize_done_
	initialize_done_ = true;

	// return value
	std::map<std::string, size_t> size_map;

	// specified values
	size_map["n_fixed"]                    = n_fixed_;
	size_map["n_random"]                   = n_random_;
	size_map["quasi_fixed"]                = size_t(  quasi_fixed_ );
	size_map["A_nr"]                       = A_rcv_.nr();
	size_map["A_nnz"]                      = A_rcv_.nnz();
	size_map["ran_like_fun.size_var"]      = ran_like_fun_.size_var();
	size_map["fix_like_fun.size_var"]      = fix_like_fun_.size_var();

	// unspecified values
	size_map["num_bytes_before"]           = num_bytes_before;
	size_map["ran_like_a1fun_.size_var()"] = ran_like_a1fun_.size_var();
	size_map["ran_hes_uu_rcv_.nnz()"]      = ran_hes_uu_rcv_.nnz();
	size_map["ran_hes_fun_.size_var()"]    = ran_hes_fun_.size_var();
	size_map["hes_cross_.subset.nnz()"]    = hes_cross_.subset.nnz();
	size_map["fix_like_jac_.subset.nnz()"] = fix_like_jac_.subset.nnz();
	size_map["fix_like_hes_.subset.nnz()"] = fix_like_hes_.subset.nnz();
	size_map["fix_con_fun_.size_var()"]    = fix_con_fun_.size_var();
	size_map["fix_con_jac_.subset.nnz()"]  = fix_con_jac_.subset.nnz();
	size_map["fix_con_hes_.subset.nnz()"]  = fix_con_hes_.subset.nnz();
	size_map["num_bytes_after"]            = CppAD::thread_alloc::inuse(thread);
	//
	if( trace_init_ )
		std::cout << "End cppad_mixed::initialize\n";
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
	catch(const std::exception& e)
	{	std::string error_message = "initialize: std::exception: ";
		error_message += e.what();
		fatal_error(error_message);
		assert(false);
	}
	catch(const CppAD::mixed::exception& e)
	{	std::string error_message = e.message("initialize");
		fatal_error(error_message);
		assert(false);
	}
# endif
	return ret;
}
