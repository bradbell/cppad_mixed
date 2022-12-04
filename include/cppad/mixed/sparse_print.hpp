// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSE_PRINT_HPP
# define CPPAD_MIXED_SPARSE_PRINT_HPP
/*
{xrst_begin sparse_print}
{xrst_spell
   cout
}

Print and Eigen Sparse Matrix
#############################

Syntax
******
``sparse_print`` ( *label* , *mat* )

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

Scalar
******
It is the type for a scalar in the matrix.
If *s* has type *Scalar* ,

   ``std::cout <<`` *s*

must print the value of *s* on standard output.

label
*****
This argument has prototype

   ``const std::string&`` *label*

Is a label printed before the matrix,
If it is empty, no label is printed.

mat
***
Is the sparse matrix which must have one of the following prototypes:

| |tab| |tab| ``const Eigen::SparseMatrix<`` *Scalar* , ``Eigen::ColMajor>&`` *mat*
| |tab| |tab| ``const Eigen::SparseMatrix<`` *Scalar* , ``Eigen::RowMajor>&`` *mat*

{xrst_end sparse_print}
*/
# include <Eigen/SparseCore>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

template <class Scalar> void sparse_print(
   const std::string&                                  label ,
   const Eigen::SparseMatrix<Scalar, Eigen::ColMajor>& mat   )
{  if( label != "" )
      std::cout << label << ":\n";
   typedef typename
   Eigen::SparseMatrix<Scalar, Eigen::ColMajor>::InnerIterator iterator;
   for(int j = 0; j < mat.outerSize(); j++)
   {  bool first = true;
      std::cout << "col " << j << ": ";
      for(iterator itr(mat,j); itr; ++itr)
      {  int i = itr.row();
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
{  if( label != "" )
      std::cout << label << ":\n";
   typedef typename
   Eigen::SparseMatrix<Scalar, Eigen::RowMajor>::InnerIterator iterator;
   for(int i = 0; i < mat.outerSize(); i++)
   {  bool first = true;
      std::cout << "rwo " << i << ": ";
      for(iterator itr(mat,i); itr; ++itr)
      {  int j = itr.col();
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
