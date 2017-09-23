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
# include <cppad/mixed/exception.hpp>

/*
$begin fix_con_jac$$
$spell
	CppAD
	cppad
	eval
	vec
	const
	Cpp
	Jacobian
	jac
$$

$section Jacobian of Fixed Constraint$$

$head Syntax$$
$icode%mixed_object%.fix_con_jac(
	%fixed_vec%, %row_out%, %col_out%, %val_out%
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
vector $latex \theta$$ at which $latex c_\theta ( \theta )$$ is evaluated.

$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to $code constraint_jac$$.
If it's input size is zero,
upon return it contains the row indices for the Jacobian elements
that are possibly non-zero.

$head col_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %col_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same as for
a previous call to $code constraint_jac$$.
If it's input size is zero,
upon return it contains the column indices for the Jacobian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$head val_out$$
This argument has prototype
$codei%
	CppAD::vector<double>& %val_out%
%$$
If the input size of this array is non-zero, it must have the same size
as for a previous call to $code constraint_jac$$.
Upon return, it contains the value of the Jacobian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$children%
	example/private/fix_con_jac.cpp
%$$
$head Example$$
The file $cref fix_con_jac.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


void cppad_mixed::fix_con_jac(
	const d_vector&        fixed_vec   ,
	CppAD::vector<size_t>& row_out     ,
	CppAD::vector<size_t>& col_out     ,
	d_vector&              val_out     )
{
	// make sure initilialize has been called
	if( ! initialize_done_  )
	{	std::string error_message =
		"fix_con_jac: initialize was not called before constraint_jac";
		fatal_error(error_message);
	}
	size_t nnz = fix_con_jac_.subset.nnz();
	if( nnz == 0 )
	{	// sparse Jacobian has no entries
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( val_out.size() == 0 );
		return;
	}
	if( row_out.size() == 0 )
	{	assert( col_out.size() == 0 );
		row_out = fix_con_jac_.subset.row();
		col_out = fix_con_jac_.subset.col();
		val_out.resize( nnz );
	}
# ifndef NDEBUG
	else
	{	for(size_t k = 0; k < nnz; k++)
		{	assert( row_out[k] == fix_con_jac_.subset.row()[k] );
			assert( col_out[k] == fix_con_jac_.subset.col()[k] );
		}
	}
# endif
	assert( row_out.size() == nnz );
	assert( col_out.size() == nnz );
	assert( val_out.size() == nnz );
	//
	sparse_rc   not_used_pattern;
	std::string not_used_coloring;
	assert( fix_con_jac_.forward == false );
	fix_con_fun_.sparse_jac_rev(
		fixed_vec            ,
		fix_con_jac_.subset  ,
		not_used_pattern     ,
		not_used_coloring    ,
		fix_con_jac_.work
	);
	for(size_t k = 0; k < nnz; k++)
		val_out[k] = fix_con_jac_.subset.val()[k];
	//
	if( CppAD::hasnan( val_out ) ) throw CppAD::mixed::exception(
		"fix_con_jac", "result has a nan"
	);
	return;
}


