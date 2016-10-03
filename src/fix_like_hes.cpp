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
$begin fix_like_hes$$
$spell
	CppAD
	cppad
	eval
	vec
	const
	Cpp
	hes
$$

$section Hessian of Fixed Likelihood$$

$head Syntax$$
$icode%mixed_object%.fix_like_hes(
	%fixed_vec%, %weight%, %row_out%, %col_out%, %val_out%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fix_likelihood$$
In the special case where $cref fix_likelihood$$ returns the empty vector,
vectors $icode row_out$$, $icode col_out$$, and $icode val_out$$ are empty.

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
components of the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
corresponding to the $cref fix_likelihood$$.
It has the same size as the corresponding return value
$cref/vec/fix_likelihood/vec/$$.

$head Hessian$$
We use $latex w$$ to denote the vector corresponding to $icode weight$$
and $latex v( \theta )$$ to denote the function corresponding th
the negative log-density vector.
The Hessian is for the function
$latex \[
	\sum_{i} w_i v_i ( \theta )
\] $$.


$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
If the input size of this array is non-zero,
the entire vector must be the same
as for a previous call to $code fix_like_hes$$.
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
a previous call to $code fix_like_hes$$.
If it's input size is zero,
upon return it contains the column indices for the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$head val_out$$
This argument has prototype
$codei%
	CppAD::vector<double>& %val_out%
%$$
If the input size of this array is non-zero, it must have the same size
as for a previous call to $code fix_like_hes$$.
Upon return, it contains the value of the Hessian elements
that are possibly non-zero (and will have the same size as $icode row_out$$).

$children%
	example/private/fix_like_hes.cpp
%$$
$head Example$$
The file $cref fix_like_hes.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


void cppad_mixed::fix_like_hes(
	const d_vector&        fixed_vec   ,
	const d_vector&        weight      ,
	CppAD::vector<size_t>& row_out     ,
	CppAD::vector<size_t>& col_out     ,
	d_vector&              val_out     )
{	assert( init_fix_like_done_ );
	assert( row_out.size() == col_out.size() );
	assert( row_out.size() == val_out.size() );
	//
	if( fix_like_hes_.row.size() == 0 )
	{	assert( fix_like_hes_.col.size() == 0 );
		assert( row_out.size() == 0 );
		assert( col_out.size() == 0 );
		assert( val_out.size() == 0 );
		return;
	}
	if( row_out.size() == 0 )
	{	row_out = fix_like_hes_.row;
		col_out = fix_like_hes_.col;
		val_out.resize( row_out.size() );
	}
# ifndef NDEBUG
	else
	{	size_t n_nonzero = fix_like_hes_.row.size();
		assert( row_out.size() == n_nonzero );
		for(size_t k = 0; k < n_nonzero; k++)
		{	assert( row_out[k] == fix_like_hes_.row[k] );
			assert( col_out[k] == fix_like_hes_.col[k] );
		}
	}
# endif

	CppAD::vector< std::set<size_t> > not_used;
	fix_like_fun_.SparseHessian(
		fixed_vec       ,
		weight          ,
		not_used        ,
		row_out         ,
		col_out         ,
		val_out         ,
		fix_like_hes_.work
	);

	return;
}


