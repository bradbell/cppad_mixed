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
$begin ranobj_hes$$
$spell
	CppAD
	ranobj
	cppad
	obj
	hes
	vec
	const
	Cpp
	xam
$$

$section Hessian of Random Part of Objective w.r.t Fixed Effects$$

$head Syntax$$
$icode%mixed_object%.ranobj_hes(
	%fixed_vec%, %random_vec%, %row_out%, %col_out%, %val_out%
)%$$

$head Purpose$$
This routine computes the Hessian w.r.t the fixed effects of
random part of the of the objective; i.e.,
$latex \[
	H_\beta^{(2)} ( \beta , \theta , u )
\] $$
see
$cref/H(beta, theta, u)
	/theory
	/Hessian of Random Objective
	/Approximate Random Objective, H(beta, theta, u)
/$$.

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
$cref/fixed effects/cppad_mixed/Fixed Effects, theta/$$
vector $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
Given the fixed effects $latex \theta$$, it is the corresponding
$cref/optimal random effects
	/theory
	/Objective
	/Optimal Random Effects, u^(theta)
/$$
$latex \hat{u} ( \theta )$$.

$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to $code ranobj_hes$$.
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
a previous call to $code ranobj_hes$$.
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
as for a previous call to $code ranobj_hes$$.
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$children%
	example/private/ranobj_hes_xam.cpp
%$$
$head Example$$
The file $cref ranobj_hes_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


// ----------------------------------------------------------------------------
// ranobj_hes
void cppad_mixed::ranobj_hes(
	const d_vector&          fixed_vec   ,
	const d_vector&          random_vec  ,
	CppAD::vector<size_t>&   row_out     ,
	CppAD::vector<size_t>&   col_out     ,
	d_vector&                val_out     )
{	assert( init_hes_ranobj_done_ );
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_vec.size() );

	// size of outputs
	size_t n_nonzero = hes_ranobj_.row.size();
	if( n_nonzero == 0 )
	{	// special case where Hessian is zero.
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( val_out.size() == 0 );
		return;
	}
	// check recording
	assert( hes_ranobj_.col.size() == n_nonzero );

	// make sure outputs have proper dimension
	assert( row_out.size() == col_out.size() );
	assert( row_out.size() == val_out.size() );

	// check if this is first call
	if( row_out.size() == 0 )
	{	row_out.resize(n_nonzero);
		col_out.resize(n_nonzero);
		val_out.resize(n_nonzero);
		for(size_t k = 0; k < n_nonzero; k++)
		{	row_out[k] = hes_ranobj_.row[k];
			col_out[k] = hes_ranobj_.col[k];
		}
	}

	// create a d_vector containing (beta, theta, u)
	d_vector beta_theta_u( 2 * n_fixed_ + n_random_ );
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);

	// create an a2d weighting vector
	d_vector w(1);
	w[0] = 1.0;

	// First call to SparseHessian is during init_hes_ranobj
	CppAD::vector< std::set<size_t> > not_used(0);

	// compute the sparse Hessian
	ranobj_fun_.SparseHessian(
		beta_theta_u,
		w,
		not_used,
		hes_ranobj_.row,
		hes_ranobj_.col,
		val_out,
		hes_ranobj_.work
	);

# ifndef NDEBUG
	for(size_t k = 0; k < n_nonzero; k++)
	{	assert( row_out[k] == hes_ranobj_.row[k] );
		assert( col_out[k] == hes_ranobj_.col[k] );
	}
# endif
}


