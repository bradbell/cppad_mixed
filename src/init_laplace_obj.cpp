// $Id:$
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
$begin init_laplace_obj$$
$spell
	init
	obj
	vec
	cppad
	hes
$$

$section Initialize Second Order of Approximate Objective and It's Hessian$$

$head Syntax$$
$icode%mixed_object%.init_laplace_obj(%fixed_vec%, %random_opt%)%$$

$head Prototype$$
$srcthisfile%
	0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1
%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head Assumptions$$
The member variable
$cref/init_ran_like_done_/init_ran_like/init_ran_like_done_/$$
is true.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which the initialization is done.

$head random_opt$$
This specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ optimization
at which the initialization is done.
It should be the optimal value given the fixed effects
so that the Hessian w.r.t the random effects is more likely to be
positive definite.

$head init_laplace_obj_done_$$
The input value of this member variable must be false.
Upon return it is true.

$head laplace_obj_fun_, init_laplace_obj_fun_done_$$
These have the same specification as in $cref init_laplace_obj_fun$$.
This initialization is done at the $icode fixed_vec$$ value for the
fixed effects and the corresponding optimal value for the random effects.
(This Hessian of the random likelihood w.r.t. the fixed effects is
more likelihood to be positive definite at the optimal random effects.)

$head laplace_obj_hes_, init_laplace_obj_hes_done_$$
These have the same specification as in $cref init_laplace_obj_hes$$.
This initialization is done at the $icode fixed_vec$$ value for the
fixed effects and the corresponding optimal value for the random effects.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

// BEGIN_PROTOTYPE
void cppad_mixed::init_laplace_obj(
	const d_vector&	   fixed_vec        ,
	const d_vector&	   random_opt       )
// END_PROTOTYPE
{	// This initialization is not called by initialize
	// so it has its own tracing.
	if( trace_init_ )
		std::cout << "Begin cppad_mixed::init_laplace_obj\n";
	//
	assert( ! init_laplace_obj_fun_done_ );
	assert( ! init_laplace_obj_hes_done_ );
	assert( ! init_laplace_obj_done_ );

	// laplace_obj_fun_
	init_laplace_obj_fun(fixed_vec, random_opt);
	assert( init_laplace_obj_fun_done_ );
	if( trace_init_ )
		std::cout << "init_laplace_obj_fun_done_\n";

	// laplace_obj_hes_
	init_laplace_obj_hes(fixed_vec, random_opt);
	assert( init_laplace_obj_hes_done_ );
	if( trace_init_ )
		std::cout << "init_laplace_obj_hes_done_\n";

	// init_laplace_obj_done_
	init_laplace_obj_done_ = true;
	if( trace_init_ )
		std::cout << "End cppad_mixed::init_laplace_obj\n";
}
