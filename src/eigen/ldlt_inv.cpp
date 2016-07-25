// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ldlt_eigen_inv$$

$section Compute a Subset of the Inverse of Factored Matrix$$
$spell
	ldlt_obj
	inv
	CppAD
	const
	eigen
$$

$head Under Construction$$

$head Syntax$$
$icode%ldlt_obj%.inv(%row_in%, %col_in%, %val_out%)%$$

$head Prototype$$
$srcfile%src/eigen/ldlt_inv.cpp
	%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
The $cref ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function solves for a subset of the inverse of the
sparse symmetric matrix that has been factored.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_eigen %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_eigen_update$$.

$head row_in$$
This vector contains the row indices for the components of the
inverse that we are computing.

$head col_in$$
This vector contains the column indices for the components of the
inverse that we are computing. It must have the same size as $icode row_in$$.

$head val_out$$
This matrix must have the same size as $icode row_in$$.
The input values of its elements do not matter.
Upon return, it
contains the values for the components of the inverse that we are computing.
To be specific, for $icode%k% = 0 , %...%, %K%-1%$$,
$icode%val_out%[%k%]%$$
is $icode%row_in%[%k%]%$$, $icode%col_in%[%k%]%$$ component of the inverse.

$head Method$$
This routine uses $cref ldlt_eigen_solve_H$$ to solve for the
requested components of the inverse one column at a time.

$end
------------------------------------------------------------------------------
*/
# include <cstddef>
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/sparseinv.hpp>
# include <cppad/utility/index_sort.hpp>

namespace CppAD { namespace mixed {

// BEGIN_PROTOTYPE
void ldlt_eigen::inv(
	const CppAD::vector<size_t>& row_in    ,
	const CppAD::vector<size_t>& col_in    ,
	CppAD::vector<double>&       val_out   )
// END_PROTOTYPE
{	using CppAD::vector;
	//
	// starting index
	size_t K  = row_in.size();
	size_t k  = 0;
	size_t row = n_row_;
	size_t col = n_row_;
	if( k < K )
	{	row = row_in[k];
		col = col_in[k];
		assert( row < n_row_ );
		assert( col < n_row_ );
	}
	//
	CppAD::vector<size_t> row_solve;
	CppAD::vector<double> rhs_solve, val_solve;
	for(size_t j = 0; j < n_row_; j++)
	{	// vectors for this column
		row_solve.resize(0);
		rhs_solve.resize(0);

		// only need for rows where f_{u,u} (theta_u) is possibly not zero
		size_t k_start        = K;
		bool   found_diagonal = false;
		while( col <= j )
		{	if( col == j )
			{	// this row of f_{u,u} (theta, u) is possibly non-zero
				row_solve.push_back( row );
				// rhs_solve needs to be j-th column of identity matrix
				// so solution is j-th column of inverse
				if( row == j )
				{	rhs_solve.push_back(1.0);
					found_diagonal = true;
				}
				else
					rhs_solve.push_back(0.0);
				//
				// index in row_in and col_in where j-th column starts
				if( k_start == K )
					k_start = k;
			}
			k++;
			if( k < K )
			{	row = row_in[k];
				col = col_in[k];
				assert( row < n_row_ );
				assert( col < n_row_ );
			}
			else
				row = col = n_row_;
		}
		assert( col > j );
		//
		// Cannot compute cholesky factor if f_{u,u} (theta, u) is zero
		if( ! found_diagonal )
		{	row_solve.push_back(j);
			rhs_solve.push_back(1.0);
		}
		// if k_start == K, we do not need any components of the inverse
		if( k_start < K )
		{	val_solve.resize( row_solve.size() );
			//
			solve_H(row_solve, rhs_solve, val_solve);
			//
			size_t nr = row_solve.size();
			if( ! found_diagonal )
				nr--;
			for(size_t ell = 0; ell < nr; ell++)
			{	size_t k    = k_start + ell;
				val_out[k] = val_solve[ell];
			}
		}
	}
	return;
}

} }
