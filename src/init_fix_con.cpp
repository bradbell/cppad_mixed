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
does not matter.
If $icode quasi_fixed$$ is false,
upon return $code fix_con_jac_$$ contains
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
does not matter.
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
// fix_con_use_set
// --------------------------------------------------------------------------
void fix_con_jac_use_set(
	const CppAD::vector<double>&       fixed_vec ,
	CppAD::ADFun<double>&              fun       ,
	CppAD::mixed::sparse_jac_info&     jac_info  )
{	jac_info.work.clear();
	assert( jac_info.row.size() == 0 );
	assert( jac_info.col.size() == 0 );
	// -----------------------------------------------------------------------
	typedef CppAD::vector< std::set<size_t> > set_sparsity;
	size_t  m = fun.Range();
	size_t  n = fun.Domain();
	assert( fixed_vec.size() == n );
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
	jac_info.direction = CppAD::mixed::sparse_jac_info::Reverse;
	// -----------------------------------------------------------------------
	// compute the work vector for reuse during future calls to
	// SparseJacobianForward
	CppAD::vector<double> not_used( jac_info.row.size() );
	fun.SparseJacobianReverse(
		fixed_vec,
		pattern,
		jac_info.row,
		jac_info.col,
		not_used,
		jac_info.work
	);
	return;
}
// --------------------------------------------------------------------------
// fix_con_use_bool
// --------------------------------------------------------------------------
void fix_con_jac_use_bool(
	const CppAD::vector<double>&       fixed_vec ,
	CppAD::ADFun<double>&              fun       ,
	CppAD::mixed::sparse_jac_info&     jac_info  )
{	jac_info.work.clear();
	assert( jac_info.row.size() == 0 );
	assert( jac_info.col.size() == 0 );
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
	jac_info.direction = CppAD::mixed::sparse_jac_info::Reverse;
	// -----------------------------------------------------------------------
	// compute the work vector for reuse during future calls to
	// SparseJacobianForward
	CppAD::vector<double> not_used( jac_info.row.size() );
	fun.SparseJacobianReverse(
		fixed_vec,
		pattern,
		jac_info.row,
		jac_info.col,
		not_used,
		jac_info.work
	);
	return;
}
} // END_EMPTY_NAMESPACE

void cppad_mixed::init_fix_con(bool bool_sparsity, const d_vector& fixed_vec  )
{	assert( fixed_vec.size() == n_fixed_ );
	assert( ! init_fix_con_done_ );

	// ------------------------------------------------------------------------
	// fix_con_fun_
	// ------------------------------------------------------------------------
	// convert to an a1d_vector
	a1d_vector a1_theta(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta[j] = fixed_vec[j];

	// start recording a1_double operations
	Independent(a1_theta);

	// compute constraint
	a1d_vector a1_vec = fix_constraint(a1_theta);
	if( a1_vec.size() == 0 )
	{	CppAD::AD<double>::abort_recording();
		init_fix_con_done_ = true;
		assert( fix_con_fun_.size_var() == 0 );
		return;
	}

	// save the recording
	fix_con_fun_.Dependent(a1_theta, a1_vec);

	// optimize the recording
# ifdef NDEBUG
	fix_con_fun_.optimize();
# endif

	// ------------------------------------------------------------------------
	// fix_con_jac_.row, fix_con_jac_.col, fix_con_jac_.work
	// ------------------------------------------------------------------------
	if( bool_sparsity) fix_con_jac_use_bool(
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
	// fix_con_hes_.row, fix_con_hes_.col, fix_con_hes_.work
	// ------------------------------------------------------------------------
	// sparsity pattern for the Hessian
	using CppAD::vectorBool;
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	size_t m            = fix_con_fun_.Range();
	size_t n            = fix_con_fun_.Domain();
	size_t n_col = vectorBool::bit_per_unit();
	sparsity_pattern pattern(n);
	vectorBool r(n * n_col), h(n * n_col);
	vectorBool s(m);
	for(size_t i = 0; i < m; i++ )
		s[i] = true;
	size_t n_loop = (n - 1) / n_col + 1;
	for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
	{	size_t j_col = i_loop * n_col;
		for(size_t i = 0; i < n; i++)
		{	for(size_t j = 0; j < n_col; j++)
				r[i * n_col + j] = (i == j_col + j);
		}
		fix_con_fun_.ForSparseJac(n_col, r);
		bool transpose = true;
		h = fix_con_fun_.RevSparseHes(n_col, s, transpose);
		//
		for(size_t i = 0; i < n; i++)
		{	for(size_t j = 0; j < n_col; j++)
			{	if( j_col + j < n )
					if( h[ i * n_col + j ] )
						pattern[i].insert(j_col + j);
			}
		}
	}
	// determine row and column indices in lower triangle of Hessian
	fix_con_hes_.row.clear();
	fix_con_hes_.col.clear();
	std::set<size_t>::iterator itr;
	for(size_t i = 0; i < n_fixed_; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			// only compute lower triangular part
			if( i >= j )
			{	fix_con_hes_.row.push_back(i);
				fix_con_hes_.col.push_back(j);
			}
		}
	}
	size_t K = fix_con_hes_.row.size();

	// compute the work vector for reuse during Hessian sparsity calculations
	d_vector weight( fix_con_fun_.Range() ), hes(K);
	fix_con_fun_.SparseHessian(
		fixed_vec       ,
		weight          ,
		pattern         ,
		fix_con_hes_.row  ,
		fix_con_hes_.col  ,
		hes             ,
		fix_con_hes_.work
	);

	init_fix_con_done_ = true;
	return;
}


