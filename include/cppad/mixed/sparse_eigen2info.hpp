// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
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
CppAD Mixed user API.

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
In this case it is assumed that there is at least one entry in every
column of $icode matrix$$.

$subhead Non-Empty on Input$$
It is assumed that on input, the size of $icode%info%.val%$$
is the same as $icode%info%.row%$$ and $icode%info%.col%$$.
Upon return, the following conditions hold:
$list number$$
The sparsity pattern in $icode info$$ is not modified.
$lnext
The value of any elements in $icode matrix$$,
that also appear in the $icode info$$, are copied to $icode info$$.
$lnext
Any elements that appear in $icode info$$, but not in $icode matrix$$,
have value zero in the return value of $icode info$$.
$lnext
The elements of $icode matrix$$,
that do not appear in $icode info$$, are ignored.
$lend

$children%example/private/sparse_eigen2info.cpp
%$$
$head Example$$
The file $cref sparse_eigen2info.cpp$$ is an example
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
			// case where input value of info is non-empty
			//
			// initilize all the values as zero
			size_t K = info.row.size();
			for(size_t k = 0; k < K; k++)
				info.val[k] = 0.0;
			//
			size_t k = 0;
			for(index j = 0; j < matrix.outerSize(); j++)
			{	// skip entries in info that are not in matrix
				size_t c = size_t(j);
				while( info.col[k] < c )
					k++;
				assert( info.col[k] == c );
				for(iterator itr(matrix, j); itr; ++itr)
				{	// skip entries in info that are not in matrix
					size_t r = size_t( itr.row() );
					while( k < K && info.col[k] == c && info.row[k] < r )
						k++;
					// check if this entry is in both matrix and info
					if( k < K && info.col[k] == c && info.row[k] == r )
						info.val[k++] = itr.value();
				}
			}
			return;
		}
} }
# endif
