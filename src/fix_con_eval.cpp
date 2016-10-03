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
$begin fix_con_eval$$
$spell
	CppAD
	cppad
	eval
	vec
	const
	Cpp
$$

$section Evaluate Fixed Constraint Function$$

$head Syntax$$
$icode%vec% = %mixed_object%.fix_con_eval(%fixed_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

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
vector $latex \theta$$ at which $latex c( \theta )$$ is evaluated.

$head vec$$
The return value has prototype
$codei%
	CppAD::vector<double> %vec%
%$$
and is the constraint function value
corresponding to the fixed effects; see
$cref/vec/fix_constraint/vec/$$.

$children%
	example/private/fix_con_eval.cpp
%$$
$head Example$$
The file $cref fix_con_eval.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


CppAD::vector<double> cppad_mixed::fix_con_eval(const d_vector& fixed_vec)
{
	// make sure initialize has been called
	if( ! initialize_done_ )
	{	std::string error_message =
		"fix_con_eval: initialize was not called before constraint_eval";
		fatal_error(error_message);
	}
	if( fix_con_fun_.size_var() == 0 )
	{	return CppAD::vector<double>(0); // empty vector
	}
	assert( fix_con_fun_.Domain() == n_fixed_ );
	return fix_con_fun_.Forward(0, fixed_vec);
}


