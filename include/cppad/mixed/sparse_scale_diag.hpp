// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_SCALE_DIAG_HPP
# define CPPAD_MIXED_SPARSE_SCALE_DIAG_HPP

/*
$begin sparse_scale_diag$$
$spell
	Eigen
	CppAD
	const
$$

$section Scales the Diagonal of an Eigen Sparse Matrix$$

$head Syntax$$
$codei%CppAD::mixed::sparse_scale_diag(%scale%, %matrix%)%$$

$head Purpose$$
Each of the diagonal elements of $icode matrix$$ is replaced by
$icode scale$$ times its original value.

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head scale$$
This argument has prototype
$codei%
	const %Float%& %scale%
%$$
where the type $icode Float$$ can be converted to the scalar
type for the matrix.

$head matrix$$
This argument has prototype
$code%
	Eigen::SparseMatrix<%Scalar%, %Options%, %Index%>& %matrix%
%$$
This is the sparse matrix for which we are scaling the diagonal elements.
The sparsity pattern for the matrix is not modified.

$subhead Scalar$$
This type must support the conversion
$codei%
	static_cast<%Scalar%>(%scale%)
%$$.

$subhead Options$$
The options for this sparse matrix are arbitrary; i.e.,
has no restrictions.

$subhead Index$$
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

$children%example/private/sparse_scale_diag.cpp
%$$
$head Example$$
The file $cref sparse_scale_diag.cpp$$ is an example
and test of $code sparse_scale_diag$$.

$end
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed {
	template <class Float, class Sparse_matrix>
	void sparse_scale_diag(Float scale, Sparse_matrix& matrix)
	{	typedef typename Sparse_matrix::Index         index;
		typedef typename Sparse_matrix::Scalar        scalar;
		typedef typename Sparse_matrix::InnerIterator iterator;
		scalar factor = static_cast<scalar>(scale);
		for(index k = 0; k < matrix.outerSize(); ++k)
		{	for(iterator itr(matrix, k); itr; ++itr)
			{	if( itr.row() == itr.col() )
					itr.valueRef() = itr.value() * scale;
			}
		}
	}
} }

# endif
