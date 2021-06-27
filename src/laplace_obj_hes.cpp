/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
$begin laplace_obj_hes$$
$spell
	rcv
	nr
	objcon
	CppAD
	ran_obj
	cppad
	obj
	hes
	vec
	const
	Cpp
	xam
$$

$section Hessian of Laplace Objective and Random Constraints$$

$head Syntax$$
$icode%mixed_object%.laplace_obj_hes(
	%fixed_vec%, %random_vec%, %weight%, %row_out%, %col_out%, %val_out%
)%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head Purpose$$
This routine computes the Hessian, w.r.t the fixed effects,
of a vector times the random part of the of the objective
and random constraints; i.e.,
$latex \[
	w_0 H_{\beta, \beta} ( \beta , \theta , u )
	+
	\sum_{i=1}^m
	w_i B^{i-1}_{\beta, \beta} ( \beta , \theta , u )
\] $$
where $latex m$$ is the number of rows in the
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$;
see
$cref/H(beta, theta, u)
	/theory
	/Approximate Laplace Objective, H(beta, theta, u)
/$$
and
$cref/B(beta, theta, u)
	/theory
	/Approximate Random Constraint Function, B(beta, theta, u)
/$$.

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
vector $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
Given the fixed effects $latex \theta$$, it is the corresponding
$cref/optimal random effects
	/theory
	/Optimal Random Effects, u^(theta)
/$$
$latex \hat{u} ( \theta )$$.

$head weight$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %weight%
%$$
It determines the weighting vector $latex w$$ in the summation
$latex \[
	w_0 H_{\beta, \beta} ( \beta , \theta , u )
	+
	\sum_{i=1}^m
	w_i B^{i-1}_{\beta, \beta} ( \beta , \theta , u )
\] $$
The size of $icode weight$$ is
$cref/A_rcv.nr()/derived_ctor/A_rcv/$$ plus one.

$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to $code laplace_obj_hes$$.
If it's input size is zero,
upon return it contains the row indices for the Hessian elements
that are possibly non-zero;
$codei%
	%row_out%[%k%] < %n_fixed%
%$$
for all $icode%k% = 0 , %...%, %row_out%.size()-1%$$

$head col_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %col_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same as for
a previous call to $code laplace_obj_hes$$.
If it's input size is zero,
upon return it contains the column indices for the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).
Note that only the lower triangle of the Hessian is computed and hence
$codei%
	%col_out%[%k%] <= %row_out%[%k%]
%$$
for all $icode%k% = 0 , %...%, %row_out%.size()-1%$$

$head val_out$$
This argument has prototype
$codei%
	CppAD::vector<double>& %val_out%
%$$
If the input size of this array is non-zero, it must have the same size
as for a previous call to $code laplace_obj_hes$$.
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$children%
	example/private/laplace_obj_hes.cpp
%$$
$head Example$$
The file $cref laplace_obj_hes.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


// ----------------------------------------------------------------------------
// laplace_obj_hes
void cppad_mixed::laplace_obj_hes(
	const d_vector&          fixed_vec   ,
	const d_vector&          random_vec  ,
	const d_vector&          weight      ,
	CppAD::vector<size_t>&   row_out     ,
	CppAD::vector<size_t>&   col_out     ,
	d_vector&                val_out     )
{	assert( init_laplace_obj_hes_done_ );
	//
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_vec.size() );
	assert( A_rcv_.nr() + 1 == weight.size() );
	//
	size_t nnz = laplace_obj_hes_.subset.nnz();
	if( nnz == 0 )
	{	// sparse Hessian has no entries
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( val_out.size() == 0 );
		return;
	}
	// check if this is first call
	if( row_out.size() == 0 )
	{	assert( col_out.size() == 0 );
		row_out = laplace_obj_hes_.subset.row();
		col_out = laplace_obj_hes_.subset.col();
		val_out.resize(nnz);
	}
# ifndef NDEBUG
	else
	{	for(size_t k = 0; k < nnz; k++)
		{	assert( row_out[k] == laplace_obj_hes_.subset.row()[k] );
			assert( col_out[k] == laplace_obj_hes_.subset.col()[k] );
		}
	}
# endif
	assert( row_out.size() == nnz );
	assert( col_out.size() == nnz );
	assert( val_out.size() == nnz );
	//
	// beta
	const d_vector& beta(fixed_vec);
	//
	// theta_u
	d_vector theta_u(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, theta_u);
	//
	// set dynamic parameters in laplace_obj_fun_
	laplace_obj_fun_.new_dynamic(theta_u);
	//
	// val_out
	sparse_rc   not_used_pattern;
	std::string not_used_coloring;
	laplace_obj_fun_.sparse_hes(
		beta                   ,
		weight                 ,
		laplace_obj_hes_.subset ,
		not_used_pattern       ,
		not_used_coloring      ,
		laplace_obj_hes_.work
	);
	for(size_t k = 0; k < nnz; k++)
		val_out[k] = laplace_obj_hes_.subset.val()[k];
	if( CppAD::hasnan( val_out ) ) CppAD::mixed::exception(
		"laplace_obj_hes", "result has a nan"
	);
	//
	return;
}
