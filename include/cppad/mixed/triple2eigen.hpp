// $Id:$
# ifndef CPPAD_MIXED_TRIPLE2EIGEN_HPP
# define CPPAD_MIXED_TRIPLE2EIGEN_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin triple2eigen$$
$spell
	Eigen
	CppAD
	nr
	nc
	const
$$

$section Convert Row, Column, Value Triple to an Eigen Sparse Matrix$$

$head Syntax$$
$icode%mat% = CppAD::mixed::triple2eigen(
	%nr%, %nc%,  %row%, %col%, %val%
)%$$

$head Private$$
This routine is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Scalar$$
Is the element type for the sparse matrix.
If $icode d$$ is a $code double$$ value,  $icode%Scalar%(d)%$$
must be a corresponding element of the matrix.

$head nr$$
This argument has prototype
$codei%
	size_t %nr%
%$$
and is the number of rows in the matrix.

$head nc$$
This argument has prototype
$codei%
	size_t %nc%
%$$
and is the number of columns in the matrix.

$head row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
and contains the row indices for possibly non-zero elements of the matrix.

$head col$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %col%
%$$
it has the same size as $icode row$$ and
contains the column indices for possibly non-zero elements of the matrix.

$head val$$
This argument has prototype
$codei%
	const CppAD::vector<%Scalar%>& %val%
%$$

$subhead Sparsity Pattern$$
If $icode%val%.size() == 0%$$,
the values in the matrix are not specified.
To be specific, for $icode%k% = 0 , %...%, %row%.size()-1%$$,
the element with index $codei%(%row%[%k%], %col%[%k%])%$$ has an unspecified
value.

$subhead Sparse Matrix$$
If $icode%val%.size() != 0%$$,
it contains the possibly non-zero values in the matrix.
To be specific, for $icode%k% = 0 , %...%, %row%.size()-1%$$,
the element with index $codei%(%row%[%k%], %col%[%k%])%$$ has value
$icode%val%[%k%]%$$.

$head mat$$
The return value has prototype
$codei%
	Eigen::SparseMatrix<%Scalar%, Eigen::ColMajor> %mat%
%$$
and is a sparse representation of the specified matrix.

$end
-----------------------------------------------------------------------------
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <class Scalar>
Eigen::SparseMatrix<Scalar> triple2eigen(
	size_t                       nr   ,
	size_t                       nc   ,
	const CppAD::vector<size_t>& row  ,
	const CppAD::vector<size_t>& col  ,
	const CppAD::vector<Scalar>& val  )
{	assert( row.size() == col.size() );
	assert( val.size() ==  0 || row.size() == val.size() );
	//
	Eigen::SparseMatrix<Scalar> ret(nr, nc);
	if( val.size() == 0 )
	{	for(size_t k = 0; k < row.size(); k++)
			ret.insert( row[k], col[k] ) = Scalar(0.0);
	}
	else
	{	for(size_t k = 0; k < row.size(); k++)
			ret.insert( row[k], col[k] ) = val[k];
	}
	//
	return ret;
}

} } // END_CPPAD_MIXED_NAMESPACE

# endif
