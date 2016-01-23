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
$begin cholmod_solve$$
$spell
	cholesky
	cholmod_obj
	const
	xam
	CppAD
$$

$section Solve Linear Equations Using Cholesky Factor$$

$head Syntax$$
$codei%%cholmod_obj%.solve(%row_in%, %val_in%, %row_out%, %val_out%)%$$

$head Private$$
The $cref cholmod$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function solve the linear equation
$latex A x = b$$ where $latex A$$ is the positive definite matrix
that has been factored,
$latex b$$ is a known column vector,
and $latex x$$ is unknown.

$head cholmod_obj$$
This object has prototype
$codei%
	const CppAD::mixed::cholmod %cholmod_obj%
%$$
In addition, it must have a previous call to
$cref cholmod_update$$.

$head row_in$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row_in%
%$$
It specifies the rows, in the column vector $latex b$$,
that are possibly non-zero.

$head val_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_in%
%$$
and it has the same size as $icode row_in$$.
It specifies the values, in the column vector $latex b$$,
that are possibly non-zero; to be specific,
for $icode%k% = 0 , %...%, %row_in%.size()-1%$$,
$codei%
	%b%[ %row_in%[%k%] ] = %val_in%[%k%]
%$$.

$head row_out$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row_out%
%$$
It specifies the rows, in the column vector $latex x$$,
that we are interested in knowing the value of.

$head val_out$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_in%
%$$
and it has the same size as $icode row_out$$.
It specifies the values, in the column vector $latex x$$,
that we are interested in.
For $icode%k% = 0 , %...%, %row_out%.size()-1%$$,
$codei%
	%b%[ %row_out%[%k%] ] = %val_out%[%k%]
%$$

$head Example$$
The file $cref/cholmod_xam.cpp/cholmod_xam.cpp/solve/$$ contains an
example and test that uses this function.

$end
*/
# include <cppad/mixed/cholmod.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cholmod::solve(
	const CppAD::vector<size_t>& row_in   ,
	const CppAD::vector<double>& val_in   ,
	const CppAD::vector<size_t>& row_out  ,
	CppAD::vector<double>&       val_out  )
{	assert( row_in.size()  == val_in.size() );
	assert( row_out.size() == val_out.size() );

	assert( rhs_ != CPPAD_NULL  );
	assert( rhs_->nrow == nrow_ );
	assert( rhs_->ncol == 1     );
	double* rhs_x = (double *) rhs_->x;
	//
# ifndef NDEBUG
	for(size_t k = 0; k < row_in.size(); k++)
		assert( row_in[k] < nrow_ );
	for(size_t k = 0; k < row_out.size(); k++)
		assert( row_out[k] < nrow_ );
	for(size_t j = 0; j < nrow_; j++)
		assert( rhs_x[j] == 0.0 );
# endif
	// set non-zero entries in right hand size
	for(size_t k = 0; k < row_in.size(); k++)
		rhs_x[ row_in[k] ] = val_in[k];

	// solve the linear equation A * sol = rhs
	int sys = CHOLMOD_A;
	// 2DO: would like to take advantate of rhs_ being sparse and only needing
	// a subset of the result vector.
	int flag = cholmod_solve2(
		sys,
		factor_,
		rhs_,
		rhs_set_,
		&sol_,
		&sol_set_,
		&work_one_,
		&work_two_,
		&common_
	);
	// check assumptions
	assert( flag == CHOLMOD_TRUE );
	assert( sol_ != CPPAD_NULL   );
	assert( sol_->nrow == nrow_ );
	assert( sol_->ncol == 1     );
	double* sol_x = (double *) sol_->x;

	// return result values
	for(size_t k = 0; k < row_out.size(); k++)
		val_out[k] = sol_x[ row_out[k] ];

	// restore the vector rhs_ to be zero
	for(size_t k = 0; k < row_in.size(); k++)
		rhs_x[ row_in[k] ] = 0.0;
}

} } // END_CPPAD_MIXED_NAMESPACE
