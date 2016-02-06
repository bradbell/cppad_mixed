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
$begin cholmod_solve2$$
$spell
	rhs
	cholesky
	cholmod_obj
	const
	xam
	CppAD
	nrow
$$

$section Solve Linear Equations Using Cholesky Factor$$

$head Syntax$$
$codei%%cholmod_obj%.solve2(%row%, %val_in%, %val_out%)%$$

$head Private$$
The $cref cholmod$$ class is an
$cref/implementation detail/cholmod/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function solves the linear equation
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

$head row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
It contains all of the rows of column vector $latex b$$ that are
non-zero and the rows of the column vector $icode x$$
that are desired.
These values must be unique; i.e.,
for $icode%j% != %k%$$,
$codei%
	%row%[%j%] != %row%[%k%]
%$$
It follows that $icode%row%.size()%$$ is less than or equal
$cref/nrow_/cholmod_ctor/nrow_/$$.

$head val_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_in%
%$$
and it has the same size as $icode row$$.
It specifies the values in the column vector $latex b$$
for each of the corresponding rows; i.e.,
for $icode%k% = 0 , %...%, %row%.size()-1%$$,
$codei%
	%b%[ %row%[%k%] ] = %val_in%[%k%]
%$$.

$head val_out$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_out%
%$$
and it has the same size as $icode row$$.
On input, the value of its elements do not matter.
Upon return, it contains the values in the column vector $latex b$$
for each of the corresponding rows; i.e.,
for $icode%k% = 0 , %...%, %row%.size()-1%$$,
$codei%
	%x%[ %row%[%k%] ] = %val_out%[%k%]
%$$.


$head Example$$
The file $cref/cholmod_xam.cpp/cholmod_xam.cpp/solve2/$$ contains an
example and test that uses this function.

$end
*/
# include <cppad/mixed/cholmod.hpp>
# include <cppad/utility/index_sort.hpp>
# include <cassert>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cholmod::solve2(
	const CppAD::vector<size_t>& row      ,
	const CppAD::vector<double>& val_in   ,
	CppAD::vector<double>&       val_out  )
{	assert( row.size() == val_in.size() );
	assert( row.size() == val_out.size() );
	assert( row.size() <= nrow_ );
# ifndef NDEBUG
	for(size_t k = 1; k < row.size(); k++)
		assert( row[k-1] < row[k] && row[k] < nrow_ );
# endif
	//
	assert( rhs_ != CPPAD_NULL  );
	assert( rhs_->nrow == nrow_ );
	assert( rhs_->ncol == 1     );
	double* rhs_x = (double *) rhs_->x;
# ifndef NDEBUG
	for(size_t j = 0; j < nrow_; j++)
		assert( rhs_x[j] == 0.0 );
# endif
	//
	assert( rhs_set_ != CPPAD_NULL );
	int* rhs_set_p = (int *) rhs_set_->p;
	int* rhs_set_i = (int *) rhs_set_->i;
	//
	// set non-zero entries in right hand size rhs_
	for(size_t k = 0; k < row.size(); k++)
		rhs_x[ row[k] ] = val_in[k];
	//
	rhs_set_p[0] = 0;
	rhs_set_p[1] = static_cast<size_t>( row.size() );
	for(size_t k = 0; k < row.size(); k++)
		rhs_set_i[k] = (int) row[k];

	// solve the linear equation A * sol = rhs
	int sys = CHOLMOD_A;
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
	//
	assert( sol_set_->nrow == nrow_ );
	assert( sol_set_->ncol == 1 );
	assert( sol_set_->xtype == CHOLMOD_PATTERN );
	assert( sol_set_->packed == CHOLMOD_TRUE);
	int* sol_set_p = (int *) sol_set_->p;
	int* sol_set_i = (int *) sol_set_->i;
	//
	assert( sol_ != CPPAD_NULL   );
	assert( sol_->nrow == nrow_ );
	assert( sol_->ncol == 1     );
	double* sol_x = (double *) sol_->x;
	//
	// sort_ is an index sort of sol_set_i
	size_t ni = (size_t) sol_set_p[1];
	key_.resize(ni);
	index_.resize(ni);
	for(size_t ell = 0; ell < ni; ell++)
		key_[ell] = (size_t) sol_set_i[ell];
	CppAD::index_sort(key_, index_);

	// return result values
	size_t k  = 0;
	for(size_t ell = 0; ell < ni; ell++)
	{	size_t i = key_[ index_[ell] ];
		assert( i <= row[k] );
		if( i == row[k] )
		{	val_out[k] = sol_x[i];
			k++;
		}
	}

	// restore the vector rhs_ to be zero
	for(size_t k = 0; k < row.size(); k++)
		rhs_x[ row[k] ] = 0.0;
}

} } // END_CPPAD_MIXED_NAMESPACE