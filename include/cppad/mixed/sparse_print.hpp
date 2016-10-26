// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_PRINT_HPP
# define CPPAD_MIXED_SPARSE_PRINT_HPP
/*
$begin sparse_print$$
$spell
	Eigen
	std
	cout
	const
$$

$section Print and Eigen Sparse Matrix$$

$head Syntax$$
$codei%sparse_print(%label%, %mat%)%$$

$head Scalar$$
Is the type for a scalar in the matrix.
If $icode s$$ has type $icode Scalar$$,
$codei%
	std::cout << %s%
%$$
must print the value of $icode s$$ on standard output.

$head label$$
Is a label printed before the matrix,
If it is empty, no label is printed.

$head mat$$
Is the sparse matrix which must have one of the following prototypes:
$codei%
		const Eigen::SparseMatrix<%Scalar%, Eigen::ColMajor>& %mat%
		const Eigen::SparseMatrix<%Scalar%, Eigen::RowMajor>& %mat%
%$$

$end
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <class Scalar> void sparse_print(
	const std::string&                                  label ,
	const Eigen::SparseMatrix<Scalar, Eigen::ColMajor>& mat   )
{	std::cout << label << ":\n";
	typedef typename
	Eigen::SparseMatrix<Scalar, Eigen::ColMajor>::InnerIterator iterator;
	for(int j = 0; j < mat.outerSize(); j++)
	{	bool first = true;
		std::cout << "col " << j << ": ";
		for(iterator itr(mat,j); itr; ++itr)
		{	int i = itr.row();
			assert( j == itr.col() );
			if( ! first )
				std::cout << ", ";
			std::cout << "(" << i << ")=" << itr.value();
			first = false;
		}
		std::cout << "\n";
	}
}

template <class Scalar> void sparse_print(
	const std::string&                                  label ,
	const Eigen::SparseMatrix<Scalar, Eigen::RowMajor>& mat   )
{	std::cout << label << ":\n";
	typedef typename
	Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator iterator;
	for(int i = 0; i < mat.outerSize(); i++)
	{	bool first = true;
		std::cout << "rwo " << i << ": ";
		for(iterator itr(mat,i); itr; ++itr)
		{	int j = itr.col();
			assert( i == itr.row() );
			if( ! first )
				std::cout << ", ";
			std::cout << "(" << j << ")=" << itr.value();
			first = false;
		}
		std::cout << "\n";
	}
}

} } // END_CPPAD_MIXED_NAMESPACE

# endif
