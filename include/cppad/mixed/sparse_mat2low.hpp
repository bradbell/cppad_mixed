// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_MAT2LOW_HPP
# define CPPAD_MIXED_SPARSE_MAT2LOW_HPP

/*
$begin sparse_mat2low$$
$spell
	Eigen
	CppAD
	const
	cols
$$

$section Extract the Lower Triangular From an Eigen Symmetric Matrix$$

$head Syntax$$
$icode%lower% = %CppAD::mixed::sparse_mat2low(%matrix%)%$$

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head matrix$$
The argument has prototype
$code%
	const Eigen::SparseMatrix<%Scalar%, %Options%, %Index%>& %matrix%
%$$
and has the same number of rows as columns; i.e.
$codei%
	%matrix%.rows() == %matrix.cols()
%$$

$subhead Scalar$$
The scalar type fro this sparse matrix are arbitrary; i.e,
has no restrictions.

$subhead Options$$
The options for this sparse matrix are arbitrary; i.e.,
has no restrictions.

$subhead Index$$
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

$head lower$$
The return value has prototype
$codei%
	Eigen::SparseMatrix<%Scalar%, %Options%, %Index%> %lower%
%$$
and is the lower triangle of $icode matrix$$; i.e.,
it has the same lower triangle as $icode matrix$$ and it has no entries
above the diagonal.

$children%example/private/sparse_mat2low_xam.cpp
%$$
$head Example$$
The file $cref sparse_mat2low_xam.cpp$$ is an example
and test of $code sparse_mat2low$$.

$end
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {

	template <class Sparse_matrix>
	Sparse_matrix sparse_mat2low(const Sparse_matrix& matrix)
	{	typedef typename Sparse_matrix::Index         index;
		typedef typename Sparse_matrix::InnerIterator iterator;
		assert( matrix.rows() == matrix.cols() );
		assert( matrix.rows() == matrix.outerSize() );
		//
		// determine the number of non-zeros entries for each outer iteration
		Eigen::Matrix<index, Eigen::Dynamic, 1> nnz( matrix.outerSize() );
		for(index k = 0; k < matrix.outerSize(); ++k)
		{	nnz[k] = index(0);
			for(iterator itr(matrix, k); itr; ++itr)
			{	if( itr.col() <= itr.row() )
					++nnz[k]; // entries is at or below the diagonal
			}
		}
		//
		// reserve space for the result
		Sparse_matrix result( matrix.rows(), matrix.rows() );
		result.reserve(nnz);
		//
		for(index k = 0; k < matrix.outerSize(); ++k)
		{	for(iterator itr(matrix, k); itr; ++itr)
			{	if( itr.col() <= itr.row() )
					result.insert(itr.row(), itr.col()) = itr.value();
			}
		}
		return result;
	}
} }

# endif
