# ifndef CPPAD_MIXED_TRIPLE2EIGEN_HPP
# define CPPAD_MIXED_TRIPLE2EIGEN_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
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
$codei%CppAD::mixed::triple2eigen(
   %mat%, %nr%, %nc%,  %row%, %col%, %val%
)%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1
%$$

$head Private$$
This routine is an implementation detail and not part of the
CppAD Mixed user API.

$head Scalar$$
is the element type for the sparse matrix.
If $icode d$$ is a $code double$$ value,  $icode%Scalar%(d)%$$
must be a corresponding element of the matrix.

$head nr$$
is the number of rows in the matrix.

$head nc$$
is the number of columns in the matrix.

$head row$$
contains the row indices for possibly non-zero elements of the matrix.

$head col$$
has the same size as $icode row$$ and
contains the column indices for possibly non-zero elements of the matrix.

$head val$$

$subhead Sparsity Pattern$$
If $icode%val%.size() == 0%$$,
the values in the matrix are not specified.
To be specific, for $icode%k% = 0 , %...%, %row%.size()-1%$$,
the element with index
$codei%(%row%[%k%], %col%[%k%])%$$ has an unspecified value.

$subhead Sparse Matrix$$
If $icode%val%.size() != 0%$$,
it contains the possibly non-zero values in the matrix.
To be specific, for $icode%k% = 0 , %...%, %row%.size()-1%$$,
the element with index $codei%(%row%[%k%], %col%[%k%])%$$ has value
$icode%val%[%k%]%$$.

$head mat$$
is a sparse representation of the specified matrix.
The input size and values in the matrix do not matter.
Upon return it is a sparse matrix with the specified size
and element values.

$end
-----------------------------------------------------------------------------
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <class Scalar>
// BEGIN_PROTOTYPE
void triple2eigen(
   Eigen::SparseMatrix<Scalar>&  mat  ,
   size_t                        nr   ,
   size_t                        nc   ,
   const s_vector&               row  ,
   const s_vector&               col  ,
   const CppAD::vector<Scalar>&  val  )
// END_PROTOTYPE
{  assert( row.size() == col.size() );
   assert( val.size() ==  0 || row.size() == val.size() );
   //
   mat.resize( int(nr), int(nc) );
   if( val.size() == 0 )
   {  for(size_t k = 0; k < row.size(); k++)
         mat.insert( int(row[k]), int(col[k]) ) = Scalar(0.0);
   }
   else
   {  for(size_t k = 0; k < row.size(); k++)
         mat.insert( int(row[k]), int(col[k]) ) = val[k];
   }
   return;
}

} } // END_CPPAD_MIXED_NAMESPACE

# endif
