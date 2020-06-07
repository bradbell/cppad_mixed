// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
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
$srcthisfile%0%// BEGIN PROTOTYPE%// END PROTOTYPE%1
%$$

$head Private$$
This function is an implementation detail and not part of the
CppAD Mixed user API.

$head left$$
This must be a square invertible lower triangular matrix; i.e.,
$icode%left%.rows() == %left%.cols()%$$ and
for each $icode%left%(%i%, %j%)%$$ in this sparse matrix,
$icode%i% >= %j%$$.

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

$children%example/private/sparse_low_tri_sol.cpp
%$$
$head Example$$
The file $cref sparse_low_tri_sol.cpp$$ is an example
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
	// One column of the result as a singly linked list with no
	// memory allocation when adding and removing elements.
	// initialize it as zero
	int res_first = nr; // null pointer value
	Eigen::Matrix<double, Dynamic, 1> res_value(nr);
	Eigen::Matrix<int, Dynamic, 1>    res_next(nr);
	res_value.setZero();
	//
	// One column of the right hand side as a singly linked list
	// with no memory allocation when adding and removing elements.
	int rhs_first;
	Eigen::Matrix<double, Dynamic, 1> rhs_value(nr);
	Eigen::Matrix<int, Dynamic, 1>    rhs_next(nr);
	//
	// for each column in the right matrix
	for(int j = 0; j < right.cols(); ++j)
	{	// iterator for non-zero entries in j-th column of right hand side
		col_itr right_itr(right, j);
		//
		// -----------------------------------------------------------------
		// convert this column of the right hand side to a singly linked list
		rhs_first = nr;   // null pointer
		if( right_itr )
		{	assert( right_itr.col() == j );
			//
			int i        = right_itr.row();
			rhs_value[i] = right_itr.value();
			rhs_first    = i;
			while( ++right_itr )
			{	assert( right_itr.col() == j );
				//
				rhs_next[i]  = right_itr.row();
				i            = right_itr.row();
				rhs_value[i] = right_itr.value();
			}
			rhs_next[i] = nr; // null pointer
		}
		// -----------------------------------------------------------------
		//
		// initialize row index for right hand side that is greater than
		// or equal the current row being solved
		int rhs_ge = rhs_first;
		//
		// initialize pointer to previous non-zero entry in result
		int res_previous = nr; // null pointer
		//
		// for each row in the left matrix
		for(int i = rhs_first; i < nr; ++i)
		{	// advance to next right hand side index greater than or equal i
			while( rhs_ge < i )
				rhs_ge = rhs_next[ rhs_ge ];
			//
			// initialize summation for (i, j) entry of the result
			if( rhs_ge == i )
			{	// --------------------------------------------------------
				// this means we have a non-zero result at index i
				if( res_previous == nr )
					res_first = i;
				else
					res_next[res_previous] = i;
				//
				res_value[i] = rhs_value[i];
				res_next[i]  = nr;
				res_previous = i;
				// --------------------------------------------------------
			}
			//
			// row rhs_first must have a non-zero result
			assert( res_previous != nr );
			//
			// for each entry in the i-th row of left matrix
			double left_ii = 0.0;
			for(row_itr left_itr(left, i); left_itr; ++left_itr)
			{	// (i, k) index in left matrix
				int k = left_itr.col();
				// check that left is lower triangular
				assert( i >= k );
				//
				if( i > k )
				{	if (res_value[k] != 0.0 )
					{	// ----------------------------------------------------
						// this means that we have a non-zero result at index i
						res_value[i] -= left_itr.value() * res_value[k];
						if( res_previous != i )
						{	res_next[res_previous] = i;
							res_next[i]            = nr;
							res_previous           = i;
						}
						// ----------------------------------------------------
					}
				}
				else
					left_ii = left_itr.value();
			}
			assert( left_ii != 0.0 );
			//
			// check if we have a non-zero entry to divide
			if( res_previous == i )
			{	res_value[i] /= left_ii;
				result.insert(i, j) = res_value[i];
			}
		}
		// restrore the res_value vector to all zeros
		int i = res_first;
		while( i < nr )
		{	res_value[i] = 0.0;
			i            = res_next[i];
		}
		// restor the res list to empty
		res_first = nr;
	}
	return result;
}
} } // END_CPPAD_MIXED_NAMESPACE