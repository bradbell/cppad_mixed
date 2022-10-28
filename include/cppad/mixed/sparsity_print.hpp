// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_SPARSITY_PRINT_HPP
# define CPPAD_MIXED_SPARSITY_PRINT_HPP
/*
$begin sparsity_print$$
$spell
	CppAD
	const
	iterator
	iterators
	std
	setvec
$$

$section Print a CppAD Internal Sparsity Pattern$$

$head Syntax$$
$codei%sparsity_print(%label%, %pattern%)%$$

$head Private$$
This routine is an implementation detail and not part of the
CppAD Mixed user API.

$head label$$
This argument has prototype
$codei%
	const std::string& %label%
%$$
It is a label printed before the matrix,
If it is empty, no label is printed.

$head pattern$$
Is the sparsity pattern which must have one of the following prototypes:
$codei%
		CppAD::local::sparse::pack_setvec& %pattern%
		CppAD::local::sparse::list_setvec& %pattern%
%$$
Note these are effectively $code const$$, but are not declared
so that the corresponding iterator can be used.
(The iterators for these types do not support $code const$$.)

$end
*/
# include <cppad/cppad.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void sparsity_print(
	const std::string&  label     ,
	CppAD::local::sparse::pack_setvec& pattern   )
{	if( label != "" )
		std::cout << label << ":\n";
	for(size_t i = 0; i < pattern.n_set(); i++)
	{	bool first = true;
			std::cout << "row " << i << ":";
		CppAD::local::sparse::pack_setvec::const_iterator itr(pattern, i);
		size_t j = *itr;
		while( j != pattern.end() )
		{	assert( j < pattern.end() );
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
{	if( label != "" )
		std::cout << label << ":\n";
	for(size_t i = 0; i < pattern.n_set(); i++)
	{	bool first = true;
			std::cout << "row " << i << ":";
		CppAD::local::sparse::list_setvec::const_iterator itr(pattern, i);
		size_t j = *itr;
		while( j != pattern.end() )
		{	assert( j < pattern.end() );
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
