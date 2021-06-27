/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_INFO2EIGEN_HPP
# define CPPAD_MIXED_SPARSE_INFO2EIGEN_HPP

/*
$begin sparse_info2eigen$$
$spell
	CppAD
	eigen
	const
	nr
	nc
$$

$section Convert a sparse_mat_info Representation to An Eigen Sparse Matrix$$

$head Syntax$$
$codei%CppAD::mixed::sparse_info2eigen(%matrix%, %info%, %nr%, %nc%)%$$

$head Private$$
This routine is an implementation detail and not part of the
CppAD Mixed user API.

$head matrix$$
This argument has prototype
$code%
	const Eigen::SparseMatrix<double, %Option%, %Index%>& %matrix%
%$$
The input value of $icode matrix$$ does not matter.
Upon return, it contains the $icode nr$$ by $icode nc$$ matrix specified by
$icode info$$.

$subhead Option$$
This must be $code Eigen::ColMajor$$.

$subhead Index$$
The index type fro this sparse matrix are arbitrary; i.e.,
has no restrictions.

$head info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %info%
%$$
This object must be in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order.

$icode nr$$
This argument has prototype
$icode%
	size_t %nr%
%$$
and is the number of rows in the matrix.

$icode nc$$
This argument has prototype
$icode%
	size_t %nc%
%$$
and is the number of columns in the matrix.

$children%example/private/sparse_info2eigen.cpp
%$$
$head Example$$
The file $cref sparse_info2eigen.cpp$$ is an example
and test of $code sparse_info2eigen$$.

$end
*/
# include <Eigen/SparseCore>
# include <cppad/mixed/sparse_mat_info.hpp>

namespace CppAD { namespace mixed {

		template <class Index>
		void sparse_info2eigen(
			Eigen::SparseMatrix<double, Eigen::ColMajor, Index>& matrix ,
			const sparse_mat_info                                info   ,
			size_t                                               nr     ,
			size_t                                               nc     )
		{	using Eigen::ColMajor;
			typedef typename
			Eigen::SparseMatrix<double, ColMajor, Index>::Index  index;
			//
			size_t total_nnz = info.col.size();
			assert( info.row.size() == total_nnz );
			//
			// compute number of non-zeros in each column
			Eigen::Matrix<index, Eigen::Dynamic, 1> nnz(nc);
			size_t k = 0;
			for(size_t j = 0; j < nc; j++)
			{	nnz[j] = 0;
				while( k < total_nnz && info.col[k] == j )
				{	nnz[j]++;
					k++;
				}
			}
			// make sure info is in column major format
			assert( k == total_nnz );
			//
			// reserve space for the matrix
			matrix.resize( int(nr), int(nc) );
			matrix.reserve(nnz);
			//
			// set the values in the matrix
			for(k = 0; k < total_nnz; k++)
			{	size_t r = info.row[k];
				size_t c = info.col[k];
				matrix.insert( int(r), int(c) ) = info.val[k];
			}
			//
			return;
		}
} }
# endif
