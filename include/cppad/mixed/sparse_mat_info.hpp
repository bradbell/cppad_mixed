// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_MAT_INFO_HPP
# define CPPAD_MIXED_SPARSE_MAT_INFO_HPP
/*
$begin sparse_mat_info$$
$spell
%$$

$section Sparse Matrix Information$$

$head Under Construction$$
This section is under construction and not yet in use.

$head Syntax$$
$codei%CppAD::mixed::sparse_mat_info %mat_info%$$

$head Purpose$$
This structure holds information about a sparse matrix.

$head row$$
The field $icode%mat_info%.row%$$ has prototype
$codei%
	CppAD::vector<size_t> %mat_info%.row
%$$
It has size zero when it is constructed.
After initialization it should contain the row indices
corresponding to possibly non-zero elements of the matrix.

$subhead K$$
We use $icode%K% = %mat_info%.row.size()%$$ below.

$head col$$
The field $icode%mat_info%.col%$$ has prototype
$codei%
	CppAD::vector<size_t> %mat_info%.col
%$$
It has size zero when it is constructed.
After initialization it should have the same size as $icode row$$
and contain the column indices
corresponding to possibly non-zero elements of the matrix.

$head val$$
The field $icode%mat_info%.val%$$ has prototype
$codei%
	CppAD::vector<double> %mat_info%.val
%$$
It has size zero when it is constructed.
After initialization it should either have size zero,
or the same size as $icode row$$
and contain the column indices
corresponding to possibly non-zero elements of the matrix.

$head Sparsity Pattern$$
If $icode%mat_info%.val.size()%$$ is zero,
we refer to $icode mat_info$$ as an sparsity pattern.
In this case, for $icode%k% = 0 , ... , %K%-1%$$,
the element with index
$codei%
	(%mat_info%.row[%k%], %mat_info%.col[%k%])
%$$
is possibly non-zero.

$head Sparse Matrix$$
If $icode%mat_info%.val.size()%$$ is equal to
$icode%mat_info%.row.size()%$$,
we refer to $icode mat_info$$ as a sparse matrix.
In this case, for $icode%k% = 0 , ... , %K%-1%$$,
the element with index
$codei%
	(%mat_info%.row[%k%], %mat_info%.col[%k%])
%$$
has value $icode%mat_info%.val[%k%]%$$ and all the
other elements of the matrix are zero.

$head Empty Matrix$$
If $icode K$$ is zero ($icode%mat_info%.row.size()%$$ is zero),
we refer to the matrix as an empty matrix.
Note that the empty matrix is the only case that is
both a sparsity pattern and a sparse matrix.

$end
*/
namespace CppAD { namespace mixed {
	struct sparse_mat_info {
		CppAD::vector<size_t>  row;
		CppAD::vector<size_t>  col;
		CppAD::vector<double>  val;
	};
} }

# endif
