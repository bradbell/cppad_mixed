// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/configure.hpp>

/*
$begin init_fix_con$$
$spell
	rcv
	Jacobians
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
$icode%mixed_object%.init_fix_con(%fixed_vec%)%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head init_fix_con_done_$$
The input value of this member variable must be false.
Upon return it is true.

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
vector $latex \theta$$ at which the initialization is done.

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
	CppAD::mixed::sparse_jac_rcv fix_con_jac_
%$$
must be empty.
Upon return, $code fix_con_jac_$$ contains the
$cref sparse_jac_rcv$$ structure for the
Jacobian of the $cref/constraints/fix_constraint/$$.

$subhead fix_con_fun_$$
This ADFun object can be used, with $code fix_con_jac_$$,
for computing sparse Jacobians; see
$cref/f/sparse_jac_rcv/Computing Sparse Jacobians/f/$$.

$head fix_con_hes_$$
The input value of
$codei%
	CppAD::mixed::sparse_hes_rcv fix_con_hes_
%$$
must be empty.
If $icode quasi_fixed$$ is false,
upon return $code fix_con_hes_$$ contains
$cref sparse_hes_rcv$$ for the
lower triangle of a weighted Hessian for the
$cref/constraints/fix_constraint/$$.

$subhead fix_con_fun_$$
This ADFun object can be used, with $code fix_con_hes_$$,
for computing sparse Hessians; see
$cref/f/sparse_hes_rcv/Computing Sparse Hessians/f/$$.

$end
*/

void cppad_mixed::init_fix_con(const d_vector& fixed_vec )
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
	size_t abort_op_index = 0;
	bool record_compare   = false;
	CppAD::Independent(a1_theta, abort_op_index, record_compare);

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
	fix_con_fun_.check_for_nan(true);

	// optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	std::string options =
		"no_conditional_skip no_compare_op no_print_for_op";
	fix_con_fun_.optimize(options);
# endif
	// ------------------------------------------------------------------------
	// fix_con_jac_
	// ------------------------------------------------------------------------
	// compute the sparsity pattern for the Jacobian
	sparse_rc pattern_in(n_fixed_, n_fixed_, n_fixed_);
	for(size_t k = 0; k < n_fixed_; k++)
		pattern_in.set(k, k, k);
	bool      transpose     = false;
	bool      dependency    = false;
	bool      internal_bool = bool_sparsity_;
	sparse_rc jac_pattern;
	fix_con_fun_.for_jac_sparsity(
		pattern_in, transpose, dependency, internal_bool, jac_pattern
	);

	// compute entire Jacobian
	fix_con_jac_.subset = d_sparse_rcv( jac_pattern );

	// use reversed mode for this sparse Jacobian
	// fix_con_jac_.Range() should be less thant fix_con_jac_.Domain()
	fix_con_jac_.forward = false;

	// use clear to make it clear that work is being computed
	// (should already be empty).
	fix_con_jac_.work.clear();

	// compute work for reuse during sparse Jacobian calculations
	std::string coloring = "cppad";
	fix_con_fun_.sparse_jac_rev(
		fixed_vec            ,
		fix_con_jac_.subset  ,
		jac_pattern          ,
		coloring             ,
		fix_con_jac_.work
	);
	// ------------------------------------------------------------------------
	// fix_con_hes_
	// ------------------------------------------------------------------------
	// no need to recalculate forward sparsity pattern.
	//
	// sparsity pattern for the Hessian
	size_t m = fix_con_fun_.Range();
	CppAD::vector<bool> select_range(m);
	for(size_t i = 0; i < m; i++)
		select_range[i] = false;
	select_range[0] = true;
	sparse_rc hes_pattern;
	fix_con_fun_.rev_hes_sparsity(
		select_range, transpose, internal_bool, hes_pattern
	);
	//
	// sparsity pattern corresponding to lower traingle of Hessian
	size_t nnz = 0;
	for(size_t k = 0; k < hes_pattern.nnz(); k++)
	{	if( hes_pattern.row()[k] >= hes_pattern.col()[k] )
			++nnz;
	}
	assert( hes_pattern.nr() == n_fixed_ );
	assert( hes_pattern.nc() == n_fixed_ );
	sparse_rc lower_pattern(n_fixed_, n_fixed_, nnz);
	size_t ell = 0;
	for(size_t k = 0; k < hes_pattern.nnz(); k++)
	{	size_t r = hes_pattern.row()[k];
		size_t c = hes_pattern.col()[k];
		if( r >= c )
			lower_pattern.set(ell++, r, c);
	}
	assert( nnz == ell );
	//
	// only compute the lower traingle of the Hessian
	fix_con_hes_.subset = d_sparse_rcv( lower_pattern );

	// compute the work vector for reuse during Hessian sparsity calculations
	// (value of weight vector does not affect work, but avoid warnings)
	d_vector weight(m);
	for(size_t i = 0; i < m; i++)
		weight[i] = 1.0;
	coloring  = "cppad.symmetric";
	fix_con_fun_.sparse_hes(
		fixed_vec            ,
		weight               ,
		fix_con_hes_.subset  ,
		hes_pattern          ,
		coloring             ,
		fix_con_hes_.work
	);
	//
	init_fix_con_done_ = true;
	return;
}
