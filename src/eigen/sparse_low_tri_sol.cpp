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
$begin sparse_low_tri_sol$$
$spell
	CppAD
	tri
	cols
$$

$section Solve a Sparse Lower Triangular Linear System$$

$head Syntax$$
$icode%result% = CppAD::mixed::sparse_low_tri_sol(%left%, %right%)%$$.

$head Prototype$$
$srcfile%src/eigen/sparse_low_tri_sol.cpp
	%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1
%$$

$head left$$
This must be a square invertible lower triangular matrix; i.e.,
$icode%left%.rows() == %left%.cols()%$$ and
for each $icode%left%(%i%, %j%)%$$ in this sparse matrix,
$icode%j% <= %i%$$.

$head right$$
The number of rows in this matrix must equal the number of columns in
$icode left$$; i.e.,
$icode%right%.rows() == %left%.cols()%$$.

$head result$$
This result has the same number of rows as $icode left$$,
the same number of columns as $icode right$$,
and satisfies the equation
$codei%
	%left% * %result% = %right%
%$$
where $code *$$ is matrix multiplication.

$children%example/private/sparse_low_tri_sol_xam.cpp
%$$
$head Example$$
The file $cref sparse_low_tri_sol_xam.cpp$$ is an example
and test of $code sparse_low_tri_sol$$.

$end
*/
# include <cppad/mixed/sparse_low_tri_sol.hpp>
# include <iostream>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

/*
template <class sparse_type>
void print(const std::string& label, const sparse_type& mat)
{	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m = mat;
	std::cout << label << "=\n" << m << "\n";
}
*/


// BEGIN PROTOTYPE
Eigen::SparseMatrix<double, Eigen::ColMajor> sparse_low_tri_sol(
	const Eigen::SparseMatrix<double, Eigen::RowMajor>&  left  ,
	const Eigen::SparseMatrix<double, Eigen::ColMajor>&  right )
// END PROTOTYPE
{	using Eigen::ColMajor;
	using Eigen::RowMajor;
	using Eigen::Dynamic;
	//
	typedef Eigen::SparseMatrix<double, ColMajor>::InnerIterator col_itr;
	typedef Eigen::SparseMatrix<double, RowMajor>::InnerIterator row_itr;
	//
	// number of rows and columns in the result
	int nr = left.rows();
	int nc = right.cols();
	assert( nr == left.cols() );
	assert( nr == right.rows() );
	//
	// initialize result matrix as empty
	Eigen::SparseMatrix<double, ColMajor> result(nr, nc);
	//
	// one column of the result as a dense vector
	Eigen::Matrix<double, Dynamic, 1> column(nr);
	//
	// for each column in the right matrix
	for(int j = 0; j < right.cols(); ++j)
	{ // iterator for non-zero entries in j-th column of right matrix
		col_itr right_itr(right, j);
		//
		// for each row in the left matrix
		for(int i = 0; i < nr; i++)
		{	//
			// initialize summation for (i, j) entry of the result
			column[i]   = 0.0;
			//
			// check if right matrix has an (i, j) entry
			if( right_itr )
			{	if( right_itr.row() == i )
				{	column[i]   = right_itr.value();
					++right_itr;
				}
			}
			//
			// for each entry in the i-th row of left matrix
			double left_ij = 0.0;
			for(row_itr left_itr(left, i); left_itr; ++left_itr)
			{	int k = left_itr.col();
				// check that left is lower triangular
				assert( k <= i );
				if( k < i )
					column[i]   -= left_itr.value() * column[k];
				else
					left_ij = left_itr.value();
			}
			assert( left_ij != 0.0 );
			// For an AD type we would need something more complicated
			// to check if the result might have been non-zero
			if( ! ( column[i] == 0.0 ) ) // check that works for nan
			{	column[i] /= left_ij;
				result.insert(i, j) = column[i];
			}
		}
	}
	return result;
}
} } // END_CPPAD_MIXED_NAMESPACE
