// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>

/*
$begin constraint_hes$$
$spell
	CppAD
	cppad
	eval
	vec
	const
	Cpp
	hes
$$

$section cppad_mixed: Hessian of Constraints w.r.t Fixed Effects$$

$head Syntax$$
$icode%mixed_object%.constraint_hes(
	%fixed_vec%, %weight%, %row_out%, %col_out%, %val_out%
)%$$

$head mixed_object$$
We use $cref/mixed_object/cppad_mixed_derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Fixed Effects, theta/$$
vector $latex \theta$$ at which the Hessian is evaluated.

$head weight$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %weight%
%$$
It specifies the value of the weights for the
components of the $cref/constraint/cppad_mixed_constraint/$$.
It has the same size as the corresponding return value
$cref/c_vec/cppad_mixed_constraint/c_vec/$$.

$head Hessian$$
We use $latex w$$ to denote the vector corresponding to $icode weight$$
and $latex c( \theta )$$ to denote the function corresponding th
the constraints
The Hessian is for the function
$latex \[
	\sum_{i} w_i c_i ( \theta )
\] $$.


$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to $code constraint_hes$$.
If it's input size is zero,
upon return it contains the row indices for the Hessian elements
that are possibly non-zero.

$head col_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %col_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same as for
a previous call to $code constraint_hes$$.
If it's input size is zero,
upon return it contains the column indices for the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$head val_out$$
This argument has prototype
$codei%
	CppAD::vector<double>& %val_out%
%$$
If the input size of this array is non-zero, it must have the same size
as for a previous call to $code constraint_hes$$.
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$children%
	example/private/constraint_hes_xam.cpp
%$$
$head Example$$
The file $cref constraint_hes_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cppad_mixed::constraint_hes(
	const d_vector&        fixed_vec   ,
	const d_vector&        weight      ,
	CppAD::vector<size_t>& row_out     ,
	CppAD::vector<size_t>& col_out     ,
	d_vector&              val_out     )
{
	// make sure initialize has been called
	if( ! initialize_done_ )
	{	std::string error_message =
		"cppad_mixed::initialize was not called before constraint_hes";
		fatal_error(error_message);
	}
	if( constraint_hes_.row.size() == 0 )
	{	// Sparse Hessian has no rows
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( constraint_hes_.row.size() == 0 );
		assert( constraint_hes_.col.size() == 0 );
		val_out.resize(0);
		return;
	}

	assert( row_out.size() == col_out.size() );
	assert( row_out.size() == val_out.size() );
	if( row_out.size() == 0 )
	{	row_out = constraint_hes_.row;
		col_out = constraint_hes_.col;
		val_out.resize( row_out.size() );
	}
# ifndef NDEBUG
	else
	{	size_t n_nonzero = constraint_hes_.row.size();
		assert( row_out.size() == n_nonzero );
		for(size_t k = 0; k < n_nonzero; k++)
		{	assert( row_out[k] == constraint_hes_.row[k] );
			assert( col_out[k] == constraint_hes_.col[k] );
		}
	}
# endif
	assert( row_out.size() != 0 );

	CppAD::vector< std::set<size_t> > not_used;
	constraint_fun_.SparseHessian(
		fixed_vec       ,
		weight          ,
		not_used        ,
		row_out         ,
		col_out         ,
		val_out         ,
		constraint_hes_.work
	);

	return;
}


} } // END_CPPAD_MIXED_NAMESPACE
