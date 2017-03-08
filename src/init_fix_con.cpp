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
# include <cppad/mixed/configure.hpp>

/*
$begin init_fix_con$$
$spell
	CppAD
	init
	cppad
	jac
	vec
	const
	Cpp
	Jacobian
	var
	hes
	bool
$$

$section Initialize Constraints as Function of Fixed Effects$$

$head Syntax$$
$icode%mixed_object%.init_fix_con(%bool_sparsity%, %fixed_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head bool_sparsity$$
This argument has prototype
$codei%
	bool %bool_sparsity%
%$$
If it is true, boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ at which the initialization is done.

$head init_fix_con_done_$$
When this function is called, this member variable must be false.
Upon return it is true.

$head fix_con_fun_$$
On input, the member variable
$codei%
	CppAD::ADFun<double> fix_con_fun_
%$$
must be empty; i.e., $code fix_con_fun_.size_var() == 0$$.
If the return value for
$cref fix_constraint$$ is empty,
$code fix_con_fun_$$ is not modified.
Otherwise,
upon return it contains the corresponding recording for the
$cref fix_constraint$$ $latex c( \theta )$$.

$head fix_con_jac_$$
The input value of
$codei%
	CppAD::mixed::sparse_jac_info fix_con_jac_
%$$
must be empty.
Upon return, $code fix_con_jac_$$ contains
$cref sparse_jac_info$$ for the
Jacobian of the $cref/constraints/fix_constraint/$$.

$subhead fix_con_fun_$$
This ADFun object can be used for the
$cref/sparse jacobian call/sparse_jac_info/Sparse Jacobian Call/f/$$.

$head fix_con_hes_$$
The input value of
$codei%
	CppAD::mixed::sparse_hes_info fix_con_hes_
%$$
must be empty.
If $icode quasi_fixed$$ is false,
upon return $code fix_con_hes_$$ contains
$cref sparse_hes_info$$ for the
lower triangle of a weighted Hessian for the
$cref/constraints/fix_constraint/$$.

$subhead fix_con_fun_$$
This ADFun object can be used for the
$cref/sparse Hessian call/sparse_hes_info/Sparse Hessian Call/f/$$.

$end
*/

namespace { // BEGIN_EMPTY_NAMESPACE
// --------------------------------------------------------------------------
// fix_con_jac_use_set
// --------------------------------------------------------------------------
void fix_con_jac_use_set(
	const CppAD::vector<double>&       fixed_vec ,
	CppAD::ADFun<double>&              fun       ,
	CppAD::mixed::sparse_jac_info&     jac_info  )
{	jac_info.work.clear();
	assert( jac_info.row.size() == 0 );
	assert( jac_info.col.size() == 0 );
	assert( jac_info.val.size() == 0 );
	// -----------------------------------------------------------------------
	typedef CppAD::vector< std::set<size_t> > set_sparsity;
	size_t  m = fun.Range();
	assert( fixed_vec.size() == fun.Domain() );
	// -----------------------------------------------------------------------
	// compute Jacobian sparsity
	// use reverse mode because m should be smaller than n
	set_sparsity r(m);
	for(size_t i = 0; i < m; i++)
		r[i].insert(i);
	set_sparsity pattern = fun.RevSparseJac(m, r);
	// -----------------------------------------------------------------------
	// Set row, column indices and direction
	std::set<size_t>::iterator itr;
	for(size_t i = 0; i < m; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			jac_info.row.push_back(i);
			jac_info.col.push_back(j);
		}
	}
	jac_info.val.resize( jac_info.row.size() );
	jac_info.direction = CppAD::mixed::sparse_jac_info::Reverse;
	// -----------------------------------------------------------------------
	// compute the work vector for reuse during future calls to
	// SparseJacobianReverse
	fun.SparseJacobianReverse(
		fixed_vec,
		pattern,
		jac_info.row,
		jac_info.col,
		jac_info.val,
		jac_info.work
	);
	return;
}
// --------------------------------------------------------------------------
// fix_con_jac_use_bool
// --------------------------------------------------------------------------
void fix_con_jac_use_bool(
	const CppAD::vector<double>&       fixed_vec ,
	CppAD::ADFun<double>&              fun       ,
	CppAD::mixed::sparse_jac_info&     jac_info  )
{	jac_info.work.clear();
	assert( jac_info.row.size() == 0 );
	assert( jac_info.col.size() == 0 );
	assert( jac_info.val.size() == 0 );
	// -----------------------------------------------------------------------
	typedef CppAD::vectorBool bool_sparsity;
	size_t  m = fun.Range();
	size_t  n = fun.Domain();
	assert( fixed_vec.size() == n );
	// -----------------------------------------------------------------------
	// compute Jacobian sparsity
	// use reverse mode because m should be smaller than n
	bool_sparsity r(m * m);
	for(size_t i = 0; i < m; i++)
	{	for(size_t j = 0; j < m; j++)
			r[i * m + j] = (i == j);
	}
	bool_sparsity pattern = fun.RevSparseJac(m, r);
	// -----------------------------------------------------------------------
	// Set row, column indices and direction
	for(size_t i = 0; i < m; i++)
	{	for(size_t j = 0; j < n; j++)
		{	if( pattern[ i * n + j ] )
			{	jac_info.row.push_back(i);
				jac_info.col.push_back(j);
			}
		}
	}
	jac_info.val.resize( jac_info.row.size() );
	jac_info.direction = CppAD::mixed::sparse_jac_info::Reverse;
	// -----------------------------------------------------------------------
	// compute the work vector for reuse during future calls to
	// SparseJacobianReverse
	fun.SparseJacobianReverse(
		fixed_vec,
		pattern,
		jac_info.row,
		jac_info.col,
		jac_info.val,
		jac_info.work
	);
	return;
}
// --------------------------------------------------------------------------
// fix_con_hes_use_set
// --------------------------------------------------------------------------
void fix_con_hes_use_set(
	const CppAD::vector<double>&       fixed_vec ,
	CppAD::ADFun<double>&              fun       ,
	CppAD::mixed::sparse_hes_info&     hes_info  )
{	hes_info.work.clear();
	assert( hes_info.row.size() == 0 );
	assert( hes_info.col.size() == 0 );
	assert( hes_info.val.size() == 0 );
	// -----------------------------------------------------------------------
	typedef CppAD::vector< std::set<size_t> > set_sparsity;
	size_t  m = fun.Range();
	size_t  n = fun.Domain();
	assert( fixed_vec.size() == n );
	// -----------------------------------------------------------------------
	// comptue Jacobian sparsity
	set_sparsity r(n);
	for(size_t i = 0; i < n; i++)
		r[i].insert(i);
	fun.ForSparseJac(n, r);
	// -----------------------------------------------------------------------
	// compute Hessian sparsity using all range components
	bool transpose = true;
	set_sparsity s(1), pattern;
	for(size_t i = 0; i < m; i++)
		s[0].insert(i);
	pattern = fun.RevSparseHes(n, s, transpose);
	// -----------------------------------------------------------------------
	// Set row, column indices
	std::set<size_t>::iterator itr;
	for(size_t i = 0; i < n; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			// only compute lower traingle
			if( j <= i )
			{	hes_info.row.push_back(i);
				hes_info.col.push_back(j);
			}
		}
	}
	hes_info.val.resize( hes_info.row.size() );
	// -----------------------------------------------------------------------
	// compute work for reuse during future calls to SparseHessian
	CppAD::vector<double> weight(m);
	for(size_t i = 0; i < m; i++)
		weight[i] = 1.0;
	fun.SparseHessian(
		fixed_vec       ,
		weight          ,
		pattern         ,
		hes_info.row    ,
		hes_info.col    ,
		hes_info.val    ,
		hes_info.work
	);
	// -----------------------------------------------------------------------
	return;
}
// --------------------------------------------------------------------------
// fix_con_hes_use_bool
// --------------------------------------------------------------------------
void fix_con_hes_use_bool(
	const CppAD::vector<double>&       fixed_vec ,
	CppAD::ADFun<double>&              fun       ,
	CppAD::mixed::sparse_hes_info&     hes_info  )
{	hes_info.work.clear();
	assert( hes_info.row.size() == 0 );
	assert( hes_info.col.size() == 0 );
	assert( hes_info.val.size() == 0 );
	// -----------------------------------------------------------------------
	typedef CppAD::vectorBool bool_sparsity;
	size_t  m = fun.Range();
	size_t  n = fun.Domain();
	assert( fixed_vec.size() == n );
	// -----------------------------------------------------------------------
	// comptue Jacobian sparsity
	bool_sparsity r(n * n);
	for(size_t i = 0; i < n; i++)
	{	for(size_t j = 0; j < n; j++)
			r[ i * n + j ] = (i == j);
	}
	fun.ForSparseJac(n, r);
	// -----------------------------------------------------------------------
	// compute Hessian sparsity using all range components
	bool transpose = true;
	bool_sparsity s(m), pattern;
	for(size_t j = 0; j < m; j++)
		s[j] = true;
	pattern = fun.RevSparseHes(n, s, transpose);
	// -----------------------------------------------------------------------
	// Set row, column indices
	for(size_t i = 0; i < n; i++)
	{	// only compute lower traingle
		for(size_t j = 0; j <= i; j++)
		{	if( pattern[i * n + j] )
			{	hes_info.row.push_back(i);
				hes_info.col.push_back(j);
			}
		}
	}
	hes_info.val.resize( hes_info.row.size() );
	// -----------------------------------------------------------------------
	// compute work for reuse during future calls to SparseHessian
	CppAD::vector<double> weight(m);
	for(size_t i = 0; i < m; i++)
		weight[i] = 1.0;
	fun.SparseHessian(
		fixed_vec       ,
		weight          ,
		pattern         ,
		hes_info.row    ,
		hes_info.col    ,
		hes_info.val    ,
		hes_info.work
	);
	// -----------------------------------------------------------------------
	return;
}
} // END_EMPTY_NAMESPACE

void cppad_mixed::init_fix_con(
	bool            bool_sparsity ,
	const d_vector& fixed_vec     )
{	assert( fixed_vec.size() == n_fixed_ );
	assert( ! init_fix_con_done_ );

	// ------------------------------------------------------------------------
	// fix_con_fun_
	// ------------------------------------------------------------------------
	// convert to an a1_vector
	a1_vector a1_theta(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta[j] = fixed_vec[j];

	// start recording a1_double operations
	Independent(a1_theta);

	// compute constraint
	a1_vector a1_vec = fix_constraint(a1_theta);
	if( a1_vec.size() == 0 )
	{	CppAD::AD<double>::abort_recording();
		init_fix_con_done_ = true;
		assert( fix_con_fun_.size_var() == 0 );
		return;
	}

	// save the recording
	fix_con_fun_.Dependent(a1_theta, a1_vec);

	// optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	fix_con_fun_.optimize();
# endif
	// ------------------------------------------------------------------------
	// fix_con_jac_
	// ------------------------------------------------------------------------
	if( bool_sparsity ) fix_con_jac_use_bool(
		fixed_vec,
		fix_con_fun_,
		fix_con_jac_
	);
	else fix_con_jac_use_set(
		fixed_vec,
		fix_con_fun_,
		fix_con_jac_
	);
	if( quasi_fixed_ )
	{	init_fix_con_done_ = true;
		return;
	}
	// ------------------------------------------------------------------------
	// fix_con_hes_
	// ------------------------------------------------------------------------
	if( bool_sparsity ) fix_con_hes_use_bool(
		fixed_vec,
		fix_con_fun_,
		fix_con_hes_
	);
	else fix_con_hes_use_set(
		fixed_vec,
		fix_con_fun_,
		fix_con_hes_
	);
	// ------------------------------------------------------------------------
	init_fix_con_done_ = true;
	return;
}


