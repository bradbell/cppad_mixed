/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
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

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
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

$head resize$$
The $code resize$$ argument has prototype
$codei%
	size_t %size%
%$$
All of the vectors,
$icode row$$, $icode col$$, and $icode val$$,
are modified to have the specified size.

$head Notation$$

$subhead Sparsity Pattern$$
We say that $icode mat_info$$ is a sparsity pattern if,
for $icode%k% = 0 , ... , %K%-1%$$,
the element with index
$codei%
	(%mat_info%.row[%k%], %mat_info%.col[%k%])
%$$
is possibly non-zero and the size or elements of
$icode%mat_info%.val%$$ are not specified.

$subhead Sparse Matrix$$
We say that $icode mat_info$$ is a sparse matrix if,
for $icode%k% = 0 , ... , %K%-1%$$,
the element with index
$codei%
	(%mat_info%.row[%k%], %mat_info%.col[%k%])
%$$
is possibly non-zero and has value $icode%mat_info%.val[%k%]%$$.

$subhead Empty Matrix$$
If $icode K$$ is zero ($icode%mat_info%.row.size()%$$ is zero),
we say that $icode mat_info$$ is the empty matrix.

$subhead Column Major Order$$
If for $icode%k% = 0 , ... , %K%-1%$$,
$codei%
	%mat_info%.col[%k%] <= %mat_info%.col[%k+1%]
	if( %mat_info%.col[%k%] == %mat_info%.col[%k+1%] )
		%mat_info%.row[%k%] < %mat_info%.row[%k+1%]
%$$
we say that $icode mat_info$$ is in column major order.

$subhead Row Major Order$$
If for $icode%k% = 0 , ... , %K%-1%$$,
$codei%
	%mat_info%.row[%k%] <= %mat_info%.row[%k+1%]
	if( %mat_info%.row[%k%] == %mat_info%.row[%k+1%] )
		%mat_info%.col[%k%] < %mat_info%.col[%k+1%]
%$$
we say that $icode mat_info$$ is in row major order.

$subhead Lower Triangular$$
If for $icode%k% = 0 , ... , %K%-1%$$,
$codei%
	%mat_info%.row[%k%] >= %mat_info%.col[%k%]
%$$
we say that $icode mat_info$$ is lower triangular.

$end
*/
# include <cppad/utility/vector.hpp>

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
