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
$begin undetermined$$
$spell
	const
	Eigen
	CppAD
	cols
	nr
	nc
$$

$section
Express An Undetermined Linear System As Dependent and Independent Variables
$$

$head Syntax$$
$icode%rank% = undetermined(%A%, %b%, %delta%, %D%, %I%, %C%, %e%)%$$

$head Purpose$$
We are give a matrix $latex A \in \B{R}^{m \times n}$$
and a vector $latex b \in \B{R}^m$$,
where $latex m < n$$ and $latex A$$ has rank $latex m$$.
Furthermore, we are interested in the linear constraint equation
$latex \[
	A x = b
\] $$
A matrix $latex C \in \B{R}^{m \times (n - m)}$$
and a vector $latex e \in \B{R}^m$$,
such that the constraint is equivalent to
$latex \[
	x_D + C x_I = e
\] $$
where $latex D$$ is a subset, of size $latex m$$,
of the column indices and $latex I$$ is the complementary subset of the
column indices.
Note that one can use this equation to
solve for the dependent  variables $latex x_D$$
as a function of independent variables $latex x_I$$.

$head A$$
This argument has prototype
$codei%
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& %A%
%$$
It is assumed that the rank of
$latex A$$ is equal to $icode%A%.rows()%$$.

$subhead nr$$
We use the notation $icode%nr% = %A%.rows()%$$; i.e.,
the number of rows in $icode A$$.

$subhead nc$$
We use the notation $icode%nc% = %A%.cols()%$$; i.e.,
the number of columns in $icode A$$.

$head b$$
This argument has prototype
$codei%
	const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& %b%
%$$
where $icode%b%.rows() = %nr%$$ and $icode%b%.cols() == 1%$$.

$head delta$$
This is the tolerance used for detecting a rank deficient matrix.
To be specific, let $latex |A|_\infty$$ be the maximum absolute element
of the matrix $latex A$$.
This algorithm uses the element with maximum absolute as a pivot for each
Gaussian elimination step.
If there is no element left with absolute value greater than
$latex \delta |A|_\infty$$,
this routine aborts with $icode%rank% !=  %nr%$$.

$head D$$
This argument has prototype
$codei%
	Eigen::Matrix<size_t, Eigen::Dynamic, 1>& %D%
%$$
where $icode%D%.rows() = %nr%$$.
The input value of its elements does not matter.
If $icode%rank% == %nr%$$,
upon return the vector $latex x_D$$ is
$codei%
	( %x%[%D%[0]] , %x%[%D%[1]] , %...% , %x%[%D%[%nr%-1]] )^T
%$$

$head I$$
This argument has prototype
$codei%
	Eigen::Matrix<size_t, Eigen::Dynamic, 1>& %I%
%$$
where $icode%I%.rows() = %nc - %nr%$$.
The input value of its elements does not matter.
If $icode%rank% == %nr%$$,
upon return the vector $latex x_I$$ is
$codei%
	( %x%[%I%[0]] , %x%[%I%[1]] , %...% , %x%[%I%[%nr%-%nc%-1]] )^T
%$$
Furthermore the union of the sets corresponding
to $latex D$$ and $latex I$$ is $latex \{ 0 , \ldots , n-1 \}$$
where $icode%n% = %nc%$$.
It follows that the sets do not intersect and none of the elements are
repeated in the vectors $icode D$$ or $icode I$$.

$head C$$
This argument has prototype
$codei%
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& %C%
%$$
and $icode%C%.rows() == %nr%$$,
$icode%C%.cols() == %nc - %nr%$$.
The input value of its elements does not matter.
If $icode%rank% == %nr%$$,
upon return it is the matrix $latex C$$ in $latex x_D + C x_I = e$$.

$head e$$
This argument has prototype
$codei%
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& %e%
%$$
and $icode%e%.rows() == %nr%$$ and $icode%e%.cols() == 1%$$.
The input value of its elements does not matter.
If $icode%rank% == %nr%$$,
upon return it is the vector $latex e$$ in $latex x_D + C x_I = e$$.

$head rank$$
The return value has prototype
$codei%
	size_t %rank%
%$$
and it is the rank of the matrix to tolerance
$cref/delta/undetermined/delta/$$.

$end
------------------------------------------------------------------------------
*/


# include <cmath>
# include <iostream>
# include <Eigen/Core>

namespace {
	using Eigen::Dynamic;
	typedef Eigen::Matrix<double, Dynamic, Dynamic> double_matrix;
	typedef Eigen::Matrix<double, Dynamic, 1>       double_vec;
	typedef Eigen::Matrix<bool,   Dynamic, 1>       bool_vec;
	typedef Eigen::Matrix<size_t, Dynamic, 1>       size_vec;
	typedef std::pair<size_t, size_t>               size_pair;
	//
	size_pair max_abs(
		const double_matrix& E        ,
		const bool_vec&      row_used ,
		const bool_vec&      col_used )
	{
		// number of rows in A
		size_t nr      = size_t( E.rows() );
		// number of columns in A
		size_t nc      = size_t( E.cols() ) - 1;
		//
		size_pair ret(0,0);
		double max_abs_val = -1.0;
		for(size_t i = 0; i < nr; i++) if( ! row_used[i] )
		{	for(size_t j = 0; j < nc; j++) if( ! col_used[j] )
			{	if( std::fabs( E(i, j) ) > max_abs_val )
				{	ret = size_pair(i, j);
					max_abs_val = std::fabs( E(i, j) );
				}
			}
		}
		// check for case where there is no possible pivot left
		assert( max_abs_val >= 0.0 );
		//
		return ret;
	}
	//
	void elementary(size_pair pivot, double_matrix& E)
	{	size_t nr = size_t( E.rows() );
		size_t r  = pivot.first;
		size_t c  = pivot.second;
		double v  = E(r, c);
		E.row(r) /= v;
		// fix roundoff on piovot element
		E(r, c) = 1.0;
		for(size_t i = 0; i < nr; i++)
		{	if( i != r )
			{	E.row(i) -= E(i, c) * E.row(r);
				// fix roundoff on pivot column
			}
		}
	}
}

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

size_t undetermined(
	const double_matrix&          A     ,
	const double_vec&             b     ,
	double                        delta ,
	size_vec&                     D     ,
	size_vec&                     I     ,
	double_matrix&                C     ,
	double_vec&                   e     )
{	size_t nr = A.rows();
	size_t nc = A.cols();
	assert(  nr < nc );
	assert(  size_t( b.rows() ) == nr );
	assert(  size_t( D.rows() ) == nr );
	assert(  size_t( I.rows() ) == nc - nr );
	assert(  size_t( C.rows() ) == nr );
	assert(  size_t( e.rows() ) == nr );
	assert(  size_t( C.cols() ) == nc - nr );
	//
	// E = [ A | b ]
	double_matrix E(nr, nc + 1 );
	E.block(0, 0,  nr, nc) = A;
	E.block(0, nc, nr, 1)  = b;
	//
	// which rows and colums have been used for pivots
	bool_vec row_used(nr), col_used(nc);
	for(size_t i = 0; i < nr; i++)
		row_used[i] = false;
	for(size_t j = 0; j < nc; j++)
		col_used[j] = false;
	//
	// mapping from pivot row to pivot column
	size_vec pivotrow2col(nr);
	//
	// maximum of absoltue value of elements of A
	double max_abs_value;
	//
	for(size_t rank = 0; rank < nr; rank++)
	{	//
		// determine the next pivot element
		size_pair pivot = max_abs(E, row_used, col_used);
		size_t r = pivot.first;
		size_t c = pivot.second;
		if( rank == 0 )
		{	// This is the element in A with maximum absolute value
			max_abs_value = std::fabs( E(r, c) );
			if( max_abs_value == 0.0 )
			{	// maxtrix has rank zero
				return rank;
			}
		}
		if( std::fabs( E(r, c) ) <= delta * max_abs_value )
			return rank;
		//
		// preform elementary row operations for this pivot
		elementary(pivot, E);
		//
		// mark this pivot row and colum as used
		row_used[r] = true;
		col_used[c] = true;
		//
		// record the column corresponding to this row
		pivotrow2col[r] = c;
	}
	//
	// D
	D = pivotrow2col;
	//
	// I and C
	size_t k = 0;
	for(size_t j = 0; j < nc; j++)
	{	// skip columns that are used for pivot operations
		if( ! col_used[j] )
		{	I[k] = j;
			C.col(k) = E.col(j);
			k++;
		}
	}
	//
	// e
	e = E.col(nc);
	//
	return nr;
}

} } // END_CPPAD_MIXED_NAMESPACE
