// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_EIGEN2INFO_HPP
# define CPPAD_MIXED_SPARSE_EIGEN2INFO_HPP

/*
$begin sparse_eigen2info$$
$spell
	CppAD
	eigen
	const
$$

$section Convert An Eigen Sparse Matrix to a sparse_mat_info Representation$$.

$head Syntax$$
$codei%CppAD::mixed::sparse_eigen2info(%matrix%, %info%)%$$

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head matrix$$
The argument has prototype
$code%
	const Eigen::SparseMatrix<double, %Option%, %Index%>& %matrix%
%$$

$subhead Option$$
This must be $code Eigen::ColMajor$$.

$subhead Index$$
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

$head info$$
This argument has prototype
$codei%
	CppAD::mixed::sparse_mat_info %info%
%$$
There are different cases depending on if it is the
$cref/empty matrix/sparse_mat_info/Notation/Empty Matrix/$$ on input;

$subhead Empty on Input$$
Upon return $icode info$$ has the same sparse matrix information
as $icode matrix$$ and is in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order.

$subhead Non-Empty on Input$$
It is assumed that on input, the size of $icode%info%.val%$$
is the same as $icode%info%.row%$$ and $icode%info%.col%$$.
Furthermore, all of the elements in $icode matrix$$ appear in the
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$ in
$icode info$$.
Upon return, this sparsity pattern is not modified.
The elements of $icode%info%.val%$$ are modified to have values corresponding
to $icode matrix$$.
Any elements that appear in $icode info$$, but not in $icode matrix$$,
have value zero in the return value of $icode info$$.

$children%example/private/sparse_eigen2info_xam.cpp
%$$
$head Example$$
The file $cref sparse_eigen2info_xam.cpp$$ is an example
and test of $code sparse_eigen2info$$.

$end
*/
# include <Eigen/SparseCore>
# include <cppad/mixed/sparse_mat_info.hpp>

namespace CppAD { namespace mixed {

		template <class Index>
		void sparse_eigen2info(
			const Eigen::SparseMatrix<double, Eigen::ColMajor, Index>& matrix ,
			sparse_mat_info&                                           info   )
		{	using Eigen::ColMajor;
			typedef typename Eigen::
			SparseMatrix<double, ColMajor, Index>::InnerIterator iterator;
			typedef typename Eigen::
			SparseMatrix<double, ColMajor, Index>::Index         index;
			//
			// case where input value of info is empty matrix
			if( info.row.size() == 0 )
			{	for(index j = 0; j < matrix.outerSize(); j++)
				{	for(iterator itr(matrix, j); itr; ++itr)
					{	info.row.push_back( size_t( itr.row() ) );
						info.col.push_back( size_t( itr.col() ) );
						info.val.push_back( itr.value() );
					}
				}
				return;
			}
			//
			// case where input value of info is non-empty
			size_t k = 0;
			for(index j = 0; j < matrix.outerSize(); j++)
			{	for(iterator itr(matrix, j); itr; ++itr)
				{	while( info.row[k] < size_t( itr.row() ) )
					{	assert( info.col[k] == size_t( itr.col() ) );
						info.val[k++] = 0.0;
					}
					assert( info.row[k] == size_t( itr.row() ) );
					assert( info.col[k] == size_t( itr.col() ) );
					info.val[k++] = itr.value();
				}
			}
			return;
		}
} }
# endif
