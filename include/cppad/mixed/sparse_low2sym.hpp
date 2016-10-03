// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_LOW2SYM_HPP
# define CPPAD_MIXED_SPARSE_LOW2SYM_HPP

/*
$begin sparse_low2sym$$
$spell
	Eigen
	CppAD
	const
	cols
	sym
$$

$section Convert an Eigen Lower Triangular Matrix To a Symmetric Matrix$$

$head Syntax$$
$icode%symmetric% = %CppAD::mixed::sparse_low2sym(%lower%)%$$

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head lower$$
The argument has prototype
$code%
	const Eigen::SparseMatrix<%Scalar%, %Options%, %Index%>& %lower%
%$$
and has the same number of rows as columns; i.e.
$codei%
	%lower%.rows() == %lower.cols()
%$$
This is a lower triangular sparse matrix; i.e.,
for each entry in $icode lower$$ has a row index
that is greater than or equal its column index.

$subhead Scalar$$
The scalar type fro this sparse matrix are arbitrary; i.e,
has no restrictions.

$subhead Options$$
The options for this sparse matrix are arbitrary; i.e.,
has no restrictions.

$subhead Index$$
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

$head symmetric$$
The return value has prototype
$codei%
	Eigen::SparseMatrix<%Scalar%, %Options%, %Index%> %symmetric%
%$$
Its lower triangle has the same entries as $icode lower$$
and it is a symmetric matrix; i.e.,
the entries above the diagonal have been set using the corresponding
entry below the diagonal.

$children%example/private/sparse_low2sym.cpp
%$$
$head Example$$
The file $cref sparse_low2sym.cpp$$ is an example
and test of $code sparse_low2sym$$.

$end
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {
	template <class Sparse_matrix>
	Sparse_matrix sparse_low2sym(const Sparse_matrix& lower)
	{	typedef typename Sparse_matrix::Index         index;
		typedef typename Sparse_matrix::InnerIterator iterator;
		assert( lower.rows() == lower.cols() );
		assert( lower.rows() == lower.outerSize() );
		//
		// determine the number of non-zeros in each row or column of result
		Eigen::Matrix<index, Eigen::Dynamic, 1> nnz( lower.outerSize() );
		for(index k = 0; k < lower.outerSize(); ++k)
		{	nnz[k] = index(0);
			for(iterator itr(lower, k); itr; ++itr)
			{	++nnz[k]; // cound entries in lower triangle
				if( itr.row() != itr.col() )
					++nnz[k]; // entries above the diagonal
			}
		}
		//
		// reserve space for the result
		Sparse_matrix result( lower.rows(), lower.rows() );
		result.reserve(nnz);
		//
		for(index k = 0; k < lower.outerSize(); ++k)
		{	for(iterator itr(lower, k); itr; ++itr)
			{	index i = itr.row();
				index j = itr.col();
				assert( j <= i );
				result.insert(i, j) = itr.value();
				if( i != j )
					result.insert(j, i) = itr.value();
			}
		}
		return result;
	}
} }

# endif
