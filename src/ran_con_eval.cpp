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
$begin ran_con_eval$$
$spell
	Au
	vec
	eval
	cppad
	const
	CppAD
$$

$section Evaluate the Random Constraint Function$$

$head Syntax$$
$icode%mixed_object%.ran_con_eval(%random_vec%, %Au%)
%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which $icode%A%*%u%$$ is evaluated.

$head Au$$
This argument has prototype
$codei%
	CppAD::vector<double> %Au%
%$$
Its size must be equal to the number of rows in the
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$.
The input value of its elements does not matter.
Upon return, it contains the product $icode%A%*%u%$$.
If the argument $icode random_vec$$ is the
$cref/optimal random effects
	/theory
	/Optimal Random Effects, u^(theta)
/$$
$icode Au$$ is the value of the
$cref/random constraint Function
	/cppad_mixed
	/Notation
	/Random Constraint Function, A*u^(theta)
/$$.

$children%
	example/private/ran_con_eval.cpp
%$$
$head Example$$
The file $cref ran_con_eval.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::ran_con_eval(
	const d_vector& random_vec ,
	d_vector&       Au         )
{	assert( random_vec.size() == n_random_ );
	assert( Au.size() == n_ran_con_ );
	assert( A_info_.row().size() == A_info_.col().size() );
	assert( A_info_.row().size() == A_info_.val().size() );
	size_t K = A_info_.row().size();
	//
	for(size_t i = 0; i < n_ran_con_; i++)
		Au[i] = 0.0;
	for(size_t k = 0; k < K; k++)
	{	size_t i = A_info_.row()[k];
		size_t j = A_info_.col()[k];
		double v = A_info_.val()[k];
		Au[i] += v * random_vec[j];
	}
	return;
}
