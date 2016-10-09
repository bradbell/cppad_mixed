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
$begin sparse_low_tri_sol.cpp$$
$spell
	tri
	cppad
$$

$section sparse_low_tri_sol: Example and Test$$

$head Private$$
This example is not part of the
$cref/cppad_mixed public API/public/$$.

$head Problem$$
This example solves the equation
$codei%
	[ 1 0 0 ]   %%         [ 1 0 0 ]
	[ 2 3 0 ] * %result% = [ 2 1 0 ]
	[ 4 5 6 ]   %%         [ 4 5 6 ]
%$$

$code
$srcfile%example/private/sparse_low_tri_sol.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$
$end
*/
// BEGIN C++
# include <cppad/mixed/sparse_low_tri_sol.hpp>

bool sparse_low_tri_sol_xam(void)
{	bool ok = true;
	double eps = 100. * std::numeric_limits<double>::epsilon();
	//
	typedef Eigen::SparseMatrix<double, Eigen::RowMajor> sparse_by_row;
	typedef Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_by_col;

	// left
	size_t nr = 3;
	sparse_by_row left(nr, nr);
	size_t count = 0;
	for(size_t i = 0; i < nr; i++)
	{	for(size_t j = 0; j <= i; j++)
			left.insert(i, j) = double( ++count );
	}
	// right
	size_t nc = 3;
	sparse_by_col right(nr, nc);
	right = left;
	//
	sparse_by_col result = CppAD::mixed::sparse_low_tri_sol(left, right);
	//
	// check that the result is a lower triangular verison of identity matrix
	for(size_t j = 0; j < nc; j++)
	{	size_t count = 0;
		for(sparse_by_col::InnerIterator itr(result, j); itr; ++itr)
		{	++count;
			double check = 0.0;
			if( itr.row() == itr.col() )
				check = 1.0;
			ok &= std::fabs( itr.value() - check ) <= eps;
			ok &= itr.col() <= itr.row();
		}
	}
	return ok;
}
// END C++