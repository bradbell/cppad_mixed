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
$begin fix_like_eval$$
$spell
	CppAD
	cppad
	eval
	vec
	const
	Cpp
$$

$section Evaluate Fixed Likelihood$$

$head Syntax$$
$icode%vec% = %mixed_object%.fix_like_eval(%fixed_vec%)%$$

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
vector $latex \theta$$ at which $latex g( \theta )$$ is evaluated.

$head vec$$
The return value has prototype
$codei%
	CppAD::vector<double> %vec%
%$$
and is a
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
corresponding to the fixed part of the negative log-likelihood
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.
To be specific;
$pre
	$$
$latex g( \theta ) = $$
$icode%vec%[0] + CppAD::abs(%vec%[1]) + %...% CppAD::abs(%vec%[%s%-1])
%$$
where $icode%s% = %vec%.size()%$$.

$children%
	example/private/fix_like_eval_xam.cpp
%$$
$head Example$$
The file $cref fix_like_eval_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


CppAD::vector<double> cppad_mixed::fix_like_eval(const d_vector& fixed_vec)
{	assert( init_fix_like_done_ );
	if( fix_likelihood_fun_.size_var() == 0 )
	{	// empty vector case
		return CppAD::vector<double>(0);
	}
	assert( fix_likelihood_fun_.Domain() == n_fixed_ );
	return fix_likelihood_fun_.Forward(0, fixed_vec);
}


