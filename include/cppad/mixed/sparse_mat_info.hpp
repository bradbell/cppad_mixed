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
	CppAD
	resize
%$$

$section Sparse Matrix Information$$

$head Syntax$$
$codei%CppAD::mixed::sparse_mat_info %mat_info
%$$
$icode%mat_info%.resize(%size%)%$$

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
or the same size as $icode row$$.

$subhead Sparsity Pattern$$
If $icode mat_info$$ is a sparsity pattern,
For $icode%k% = 0 , ... , %K%-1%$$,
the element with index
$codei%
	(%mat_info%.row[%k%], %mat_info%.col[%k%])
%$$
is possibly non-zero.
In this case, the elements of the vector $icode%mat_info%.val%$$ are
not specified.

$subhead Sparse Matrix$$
If $icode mat_info$$ is a sparse matrix,
For $icode%k% = 0 , ... , %K%-1%$$,
the element with index
$codei%
	(%mat_info%.row[%k%], %mat_info%.col[%k%])
%$$
is possibly non-zero and has value $icode%mat_info%.val[%k%]%$$.

$subhead Empty Matrix$$
If $icode K$$ is zero ($icode%mat_info%.row.size()%$$ is zero),
we refer to the matrix as the empty matrix.

$head resize$$
The $code resize$$ argument has prototype
$codei%
	size_t %size%
%$$
All of the vectors,
$code row$$, $code col$$, and $code val$$,
are modified to have the specified size.

$end
*/
namespace CppAD { namespace mixed {
	struct sparse_mat_info {
		CppAD::vector<size_t>  row;
		CppAD::vector<size_t>  col;
		CppAD::vector<double>  val;
		//
		void resize(size_t size)
		{	row.resize(size);
			col.resize(size);
			val.resize(size);
		}
	};
} }

# endif
