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
$codei%info% = CppAD::mixed::sparse_eigen2info(%matrix%)%$$

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head matrix$$
The argument has prototype
$code%
	const Eigen::SparseMatrix<double, %Option%, %Index%>& %matrix%
%$$

$subhead Option$$
This must be either $code Eigen::RowMajor$$ or $code Eigen::ColMajor$$.

$subhead Index$$
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

$head info$$
The return value has prototype
$codei%
	CppAD::mixed::sparse_mat_info %info%
%$$
It has the same sparse matrix information as $icode matrix$$ and is in
$cref/row major/sparse_mat_info/Notation/Row Major Order/$$ or
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$
depending on the value of $icode Option$$.

$head Sparsity Pattern$$
In the special case where all of the entries in $icode matrix$$
have the value $code nan$$ (not a number),
$icode%info%.val.size()%$$ is zero and $icode info$$ is a
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$.

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
		template <int Option, class Index>
		sparse_mat_info sparse_eigen2info(
			const Eigen::SparseMatrix<double, Option, Index>& matrix )
		{
			typedef typename
			Eigen::SparseMatrix<double, Option, Index>::InnerIterator iterator;
			typedef typename
			Eigen::SparseMatrix<double, Option, Index>::Index         index;
			//
			sparse_mat_info result;
			bool all_nan = true;
			for(index j = 0; j < matrix.outerSize(); j++)
			{	for(iterator itr(matrix, j); itr; ++itr)
				{	result.row.push_back( size_t( itr.row() ) );
					result.col.push_back( size_t( itr.col() ) );
					result.val.push_back( itr.value() );
					// check if itr.value() is nan
					all_nan &= itr.value() != itr.value();
				}
			}
			if( all_nan )
				result.val.resize(0);
			return result;
		}
} }
# endif
