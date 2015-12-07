// $Id:$
/* --------------------------------------------------------------------------
dismod_at: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <dismod_at/cppad_mixed.hpp>

/*
$begin init_constraint$$
$spell
	init
	cppad
	jac
	vec
	const
	Cpp
	Jacobian
	var
	hes
$$

$section Initialize Constraints as Function of Fixed Effects$$

$head Syntax$$
$codei%init_constraint(%fixed_vec%)%$$

$head Private$$
This function is $code private$$ to the $code cppad_mixed$$ class
and cannot be used by a derived
$cref/mixed_object/cppad_mixed_derived_ctor/mixed_object/$$.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Fixed Effects, theta/$$
vector $latex \theta$$ at which the initialization is done.

$head init_constraint_done_$$
When this function is called, this member variable must be false.
Upon return it is true.

$head constraint_fun_$$
On input, the member variable
$codei%
	CppAD::ADFun<double> constraint_fun_
%$$
must be empty; i.e., $code constraint_fun_.size_var() == 0$$.
If the return value for
$cref/constraint/cppad_mixed_constraint/$$ is empty,
$code constraint_fun_$$ is not modified.
Otherwise,
upon return it contains the corresponding recording for the
$cref/constraint/cppad_mixed_constraint/$$ $latex c( \theta )$$.

$head constraint_jac_$$
The input value of
$codei%
	sparse_jac_info constraint_jac_
%$$
does not matter.
If $icode quasi_fixed$$ is false,
upon return $code constraint_jac_$$ contains
$cref sparse_jac_info$$ for the
Jacobian of the $cref/constraints/cppad_mixed_constraint/$$.

$subhead constraint_fun_$$
This ADFun object can be used for the
$cref/sparse jacobian call/sparse_jac_info/Sparse Jacobian Call/f/$$.

$head constraint_hes_$$
The input value of
$codei%
	sparse_hes_info constraint_hes_
%$$
does not matter.
If $icode quasi_fixed$$ is false,
upon return $code constraint_hes_$$ contains
$cref sparse_hes_info$$ for the
lower triangle of a weighted Hessian for the
$cref/constraints/cppad_mixed_constraint/$$.

$subhead constraint_fun_$$
This ADFun object can be used for the
$cref/sparse Hessian call/sparse_hes_info/Sparse Hessian Call/f/$$.

$end
*/
namespace dismod_at { // BEGIN_DISMOD_AT_NAMESPACE

void cppad_mixed::init_constraint(const d_vector& fixed_vec  )
{	assert( fixed_vec.size() == n_fixed_ );
	assert( ! init_constraint_done_ );

	// ------------------------------------------------------------------------
	// constraint_fun_
	// ------------------------------------------------------------------------
	// convert to an a1d_vector
	a1d_vector a1_theta(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta[j] = fixed_vec[j];

	// start recording a1_double operations
	Independent(a1_theta);

	// compute constraint
	a1d_vector a1_vec = constraint(a1_theta);
	if( a1_vec.size() == 0 )
	{	CppAD::AD<double>::abort_recording();
		init_constraint_done_ = true;
		assert( constraint_fun_.size_var() == 0 );
		return;
	}

	// save the recording
	constraint_fun_.Dependent(a1_theta, a1_vec);

	// optimize the recording
# ifdef NDEBUG
	constraint_fun_.optimize();
# endif

	// ------------------------------------------------------------------------
	// constraint_jac_.row, constraint_jac_.col, constraint_jac_.work
	// ------------------------------------------------------------------------
	// compute the sparsity pattern for the Jacobian
	using CppAD::vectorBool;
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	size_t m            = constraint_fun_.Range();
	size_t n            = constraint_fun_.Domain();
	sparsity_pattern pattern(m);
	if( n < m )
	{	size_t n_col = vectorBool::bit_per_unit();
		vectorBool s(m * n_col), r(n * n_col);
		size_t n_loop = (n - 1) / n_col + 1;
		for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
		{	size_t j_col = i_loop * n_col;
			for(size_t i = 0; i < n; i++)
			{	for(size_t j = 0; j < n_col; j++)
					r[i * n_col + j] = (i == j_col + j);
			}
			s = constraint_fun_.ForSparseJac(n_col, r);
			for(size_t i = 0; i < m; i++)
			{	for(size_t j = 0; j < n_col; j++)
				{	if( j_col + j < n )
						if( s[ i * n_col + j ] )
							pattern[i].insert(j_col + j);
				}
			}
		}
	}
	else
	{	size_t n_row = vectorBool::bit_per_unit();
		vectorBool s(n_row * n), r(n_row * m);
		size_t n_loop = (n - 1) / n_row + 1;
		for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
		{	size_t i_row = i_loop * n_row;
			for(size_t i = 0; i < n_row; i++)
			{	for(size_t j = 0; j < m; j++)
					r[i * m + j] = (i_row + i == j);
			}
			s = constraint_fun_.RevSparseJac(n_row, r);
			for(size_t i = 0; i < n_row; i++)
			{	for(size_t j = 0; j < n; j++)
				{	if( i_row + i < m )
						if( s[ i * n + j ] )
							pattern[i_row + i].insert(j);
				}
			}
		}
	}
	// convert sparsity to row and column index form
	constraint_jac_.row.clear();
	constraint_jac_.col.clear();
	std::set<size_t>::iterator itr;
	for(size_t i = 0; i < m; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			constraint_jac_.row.push_back(i);
			constraint_jac_.col.push_back(j);
		}
	}

	// set direction for the sparse jacobian calculation
	constraint_jac_.direction = sparse_jac_info::Forward;

	// compute the work vector for reuse during Jacobian sparsity calculations
	d_vector jac( constraint_jac_.row.size() );
	constraint_fun_.SparseJacobianForward(
		fixed_vec       ,
		pattern         ,
		constraint_jac_.row  ,
		constraint_jac_.col  ,
		jac             ,
		constraint_jac_.work
	);
	if( quasi_fixed_ )
	{	init_constraint_done_ = true;
		return;
	}
	// ------------------------------------------------------------------------
	// constraint_hes_.row, constraint_hes_.col, constraint_hes_.work
	// ------------------------------------------------------------------------
	// sparsity pattern for the Hessian
	size_t n_col = vectorBool::bit_per_unit();
	pattern.clear();
	pattern.resize(n);
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
		constraint_fun_.ForSparseJac(n_col, r);
		bool transpose = true;
		h = constraint_fun_.RevSparseHes(n_col, s, transpose);
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
	constraint_hes_.row.clear();
	constraint_hes_.col.clear();
	for(size_t i = 0; i < n_fixed_; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			// only compute lower triangular part
			if( i >= j )
			{	constraint_hes_.row.push_back(i);
				constraint_hes_.col.push_back(j);
			}
		}
	}
	size_t K = constraint_hes_.row.size();

	// compute the work vector for reuse during Hessian sparsity calculations
	d_vector weight( constraint_fun_.Range() ), hes(K);
	constraint_fun_.SparseHessian(
		fixed_vec       ,
		weight          ,
		pattern         ,
		constraint_hes_.row  ,
		constraint_hes_.col  ,
		hes             ,
		constraint_hes_.work
	);

	init_constraint_done_ = true;
	return;
}


} // END_DISMOD_AT_NAMESPACE
