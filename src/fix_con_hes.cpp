// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>

/*
$begin fix_con_hes$$
$spell
	CppAD
	cppad
	eval
	vec
	const
	Cpp
	hes
$$

$section Hessian of Fixed Constraints$$

$head Syntax$$
$icode%mixed_object%.fix_con_hes(
	%fixed_vec%, %weight%, %row_out%, %col_out%, %val_out%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which the Hessian is evaluated.

$head weight$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %weight%
%$$
It specifies the value of the weights for the
components of the $cref fix_constraint$$.
It has the same size as the corresponding return value
$cref/vec/fix_constraint/vec/$$.

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
	example/private/fix_con_hes.cpp
%$$
$head Example$$
The file $cref fix_con_hes.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


void cppad_mixed::fix_con_hes(
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
	size_t nnz = fix_con_hes_.subset.nnz();
	if( nnz == 0 )
	{	// Sparse Hessian has no entries
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( val_out.size() == 0 );
		return;
	}
	if( row_out.size() == 0 )
	{	assert( col_out.size() == 0 );
		row_out = fix_con_hes_.subset.row();
		col_out = fix_con_hes_.subset.col();
		val_out.resize(nnz);
	}
# ifndef NDEBUG
	else
	{	for(size_t k = 0; k < nnz; k++)
		{	assert( row_out[k] == fix_con_hes_.subset.row()[k] );
			assert( col_out[k] == fix_con_hes_.subset.col()[k] );
		}
	}
# endif
	assert( row_out.size() == nnz );
	assert( col_out.size() == nnz );
	assert( val_out.size() == nnz );
	//
	sparse_rc   not_used_pattern;
	std::string not_used_coloring;
	fix_con_fun_.sparse_hes(
		fixed_vec            ,
		weight               ,
		fix_con_hes_.subset  ,
		not_used_pattern     ,
		not_used_coloring    ,
		fix_con_hes_.work
	);
	for(size_t k = 0; k < nnz; k++)
		val_out[k] = fix_con_hes_.subset.val()[k];
	//
	return;
}


