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
$begin ran_objcon_hes$$
$spell
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

$section Hessian of Random Objective and Random Constraints$$

$head Syntax$$
$icode%mixed_object%.ran_objcon_hes(
	%fixed_vec%, %random_vec%, %weight%, %row_out%, %col_out%, %val_out%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

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
	/Approximate Random Objective, H(beta, theta, u)
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
$cref/n_ran_con_/init_ran_con/n_ran_con_/$$ plus one.

$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to $code ran_objcon_hes$$.
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
a previous call to $code ran_objcon_hes$$.
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
as for a previous call to $code ran_objcon_hes$$.
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$children%
	example/private/ran_objcon_hes_xam.cpp
%$$
$head Example$$
The file $cref ran_objcon_hes_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


// ----------------------------------------------------------------------------
// ran_objcon_hes
void cppad_mixed::ran_objcon_hes(
	const d_vector&          fixed_vec   ,
	const d_vector&          random_vec  ,
	const d_vector&          weight      ,
	CppAD::vector<size_t>&   row_out     ,
	CppAD::vector<size_t>&   col_out     ,
	d_vector&                val_out     )
{	assert( init_ran_objcon_hes_done_ );
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_vec.size() );
	assert( n_ran_con_ + 1 == weight.size() );

	// size of outputs
	size_t n_nonzero = ran_objcon_hes_.row.size();
	if( n_nonzero == 0 )
	{	// special case where Hessian is zero.
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( val_out.size() == 0 );
		return;
	}
	// check recording
	assert( ran_objcon_hes_.col.size() == n_nonzero );

	// make sure outputs have proper dimension
	assert( row_out.size() == col_out.size() );
	assert( row_out.size() == val_out.size() );

	// check if this is first call
	if( row_out.size() == 0 )
	{	row_out.resize(n_nonzero);
		col_out.resize(n_nonzero);
		val_out.resize(n_nonzero);
		for(size_t k = 0; k < n_nonzero; k++)
		{	row_out[k] = ran_objcon_hes_.row[k];
			col_out[k] = ran_objcon_hes_.col[k];
		}
	}

	// create a d_vector containing (beta, theta, u)
	d_vector beta_theta_u( 2 * n_fixed_ + n_random_ );
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);

	// First call to SparseHessian is during init_ran_objcon_hes
	CppAD::vector< std::set<size_t> > not_used(0);

	// compute the sparse Hessian
	ran_objcon_fun_.SparseHessian(
		beta_theta_u,
		weight,
		not_used,
		ran_objcon_hes_.row,
		ran_objcon_hes_.col,
		val_out,
		ran_objcon_hes_.work
	);

# ifndef NDEBUG
	for(size_t k = 0; k < n_nonzero; k++)
	{	assert( row_out[k] == ran_objcon_hes_.row[k] );
		assert( col_out[k] == ran_objcon_hes_.col[k] );
	}
# endif
}
