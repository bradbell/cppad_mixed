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
$begin init_ran_con$$
$spell
	ldlt_eigen
	init
	cppad
	const
	CppAD
	cholesky
	eigen
$$

$section Initialize Random Constraints$$

$head Syntax$$
$icode%mixed_object%.init_ran_con()
%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head n_ran_con_$$
The input value of the member variable
$codei%
	size_t n_ran_con_
%$$
does not matter.
Upon return it is the number of random constraints; i.e.,
the number of rows in the matrix
$cref/A_info/derived_ctor/A_info/$$.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_ran_con(void)
{	assert( ! init_ran_con_done_ );

	// number of possibly non-zero entries
	size_t K = A_info_.row.size();
	//
	assert( n_random_ > 0 );
	assert( A_info_.col.size() == K );
	assert( A_info_.val.size() == K );

	// determine maximum row and column index
	size_t max_row_p1 = 0;
	size_t max_col_p1 = 0;
	for(size_t k = 0; k < K; k++)
	{	assert( A_info_.val[k] != 0.0 );
		max_row_p1 = std::max(max_row_p1, A_info_.row[k] + 1);
		max_col_p1 = std::max(max_col_p1,    A_info_.col[k] + 1);
	}
	assert( max_col_p1 <= n_random_ );

	// set the number of random constraints
	n_ran_con_ = max_row_p1;

	init_ran_con_done_ = true;
	return;
}
