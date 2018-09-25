// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
# include <cppad/mixed/configure.hpp>

/*
$begin init_ran_like$$
$spell
	CppAD
	init
	cppad
	vec
	const
	Cpp
$$

$section Initialize Random Likelihood$$

$head Syntax$$
$icode%mixed_object%.init_ran_like(%fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head init_ran_like_done_$$
The input value of this member variable must be false.
Upon return it is true.

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
vector $latex \theta$$ at which the initialization is done.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which the initialization is done.

$head ran_like_fun_$$
The input value of the member variable
$codei%
	CppAD::ADFun<double> ran_like_fun_
%$$
does not matter.
Upon return it contains a recording of the function
$cref ran_likelihood$$.

$head ran_like_a1fun_$$
The input value of the member variable
$codei%
	CppAD::ADFun<double> ran_like_a1fun_
%$$
does not matter.
Upon return it contains a recording of the function
$cref ran_likelihood$$.

$end
*/

void cppad_mixed::init_ran_like(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_ran_like_done_ );
	//
	using CppAD::AD;
	using CppAD::ADFun;
	using CppAD::vector;
	using CppAD::Independent;
	//
	// ------------------------------------------------------------------
	// record ran_like_a2fun
	// ------------------------------------------------------------------
	// combine into one vector
	a3_vector a3_both( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, a3_both);

	// start recording a3_double operations
	size_t abort_op_index = 0;
	bool record_compare   = false;
	Independent(a3_both, abort_op_index, record_compare);

	// extract the fixed and random effects
	a3_vector a3_theta(n_fixed_), a3_u(n_random_);
	unpack(a3_theta, a3_u, a3_both);

	// compute ran_likelihood using a3_double operations
	a3_vector a3_vec = ran_likelihood(a3_theta, a3_u);
	if( a3_vec.size() == 0 )
	{	std::string error_message =
			"init_ran_like: n_random > 0 and ran_likelihood has size 0";
		fatal_error(error_message);
	}
	if( a3_vec.size() != 1 )
	{	std::string error_message =
		"init_ran_like: ran_likelihood does not have size zero or one.";
		fatal_error(error_message);
	}

	// save the recording
	ADFun<a2_double> ran_like_a2fun(a3_both, a3_vec);
	ran_like_a2fun.check_for_nan(false);

	// optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	ran_like_a2fun.optimize("no_conditional_skip");
# endif
	// ------------------------------------------------------------------
	// record ran_like_a1fun_
	// ------------------------------------------------------------------
	// combine into one vector
	a2_vector a2_both( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, a2_both);

	// start recording a2_double operations
	Independent(a2_both, abort_op_index, record_compare);

	// compute ran_likelihood using a2_double operations
	a2_vector a2_vec = ran_like_a2fun.Forward(0, a2_both);
	assert( a2_vec.size() > 0 );
	assert( a2_vec.size() == 1 );

	// save the recording
	ran_like_a1fun_.Dependent(a2_both, a2_vec);
	ran_like_a1fun_.check_for_nan(false);

	// Brad thinks that re-optimizing will not help

	// ------------------------------------------------------------------
	// record ran_like_fun_
	// ------------------------------------------------------------------
	// combine into one vector
	a1_vector a1_both( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, a1_both);

	// start recording a1_double operations
	Independent(a1_both, abort_op_index, record_compare);

	// compute ran_likelihood using a1_double operations
	a1_vector a1_vec = ran_like_a1fun_.Forward(0, a1_both);
	assert( a1_vec.size() > 0 );
	if( a1_vec.size() == 1 );

	// save the recording
	ran_like_fun_.Dependent(a1_both, a1_vec);
	ran_like_fun_.check_for_nan(false);

	// ------------------------------------------------------------------
	init_ran_like_done_ = true;
	return;
}
