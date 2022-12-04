// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSITY_PRINT_HPP
# define CPPAD_MIXED_SPARSITY_PRINT_HPP
/*
{xrst_begin sparsity_print}
{xrst_spell
   iterator
   iterators
   setvec
}

Print a CppAD Internal Sparsity Pattern
#######################################

Syntax
******
``sparsity_print`` ( *label* , *pattern* )

Private
*******
This routine is an implementation detail and not part of the
CppAD Mixed user API.

label
*****
This argument has prototype

   ``const std::string&`` *label*

It is a label printed before the matrix,
If it is empty, no label is printed.

pattern
*******
Is the sparsity pattern which must have one of the following prototypes:

| |tab| |tab| ``CppAD::local::sparse::pack_setvec&`` *pattern*
| |tab| |tab| ``CppAD::local::sparse::list_setvec&`` *pattern*

Note these are effectively ``const`` , but are not declared
so that the corresponding iterator can be used.
(The iterators for these types do not support ``const`` .)

{xrst_end sparsity_print}
*/
# include <cppad/cppad.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void sparsity_print(
   const std::string&  label     ,
   CppAD::local::sparse::pack_setvec& pattern   )
{  if( label != "" )
      std::cout << label << ":\n";
   for(size_t i = 0; i < pattern.n_set(); i++)
   {  bool first = true;
         std::cout << "row " << i << ":";
      CppAD::local::sparse::pack_setvec::const_iterator itr(pattern, i);
      size_t j = *itr;
      while( j != pattern.end() )
      {  assert( j < pattern.end() );
         if( ! first )
            std::cout << ",";
         std::cout << " " << j;
         first = false;
         j     = *(++itr);
      }
      std::cout << "\n";
   }
}

void sparsity_print(
   const std::string&  label     ,
   CppAD::local::sparse::list_setvec& pattern   )
{  if( label != "" )
      std::cout << label << ":\n";
   for(size_t i = 0; i < pattern.n_set(); i++)
   {  bool first = true;
         std::cout << "row " << i << ":";
      CppAD::local::sparse::list_setvec::const_iterator itr(pattern, i);
      size_t j = *itr;
      while( j != pattern.end() )
      {  assert( j < pattern.end() );
         if( ! first )
            std::cout << ",";
         std::cout << " " << j;
         first = false;
         j     = *(++itr);
      }
      std::cout << "\n";
   }
}


} } // END_CPPAD_MIXED_NAMESPACE

# endif
