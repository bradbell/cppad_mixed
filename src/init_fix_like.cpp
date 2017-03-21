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
$begin init_fix_like$$
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
$$

$section Initialize Fixed Likelihood$$

$head Syntax$$
$icode%mixed_object%.init_fix_like(%fixed_vec%)%$$

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
vector $latex \theta$$ at which the initialization is done.

$head init_fix_like_done_$$
When this function is called, this member variable must be false.
Upon return it is true.

$head fix_like_fun_$$
On input, the member variable
$codei%
	CppAD::ADFun<double> fix_like_fun_
%$$
must be empty; i.e., $code fix_like_fun_.size_var() == 0$$.
If the return value for
$cref fix_likelihood$$ is empty,
$code fix_like_fun_$$ is not modified.
Otherwise,
upon return it contains the corresponding recording for the
$cref fix_likelihood$$.
The function result is the
$cref/negative log-density vector/cppad_mixed/Negative Log-Density Vector/$$
corresponding to the function
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.

$head fix_like_jac_$$
The input value of
$codei%
	CppAD::mixed::sparse_jac_rcv fix_like_jac_
%$$
must be empty.
If the return value for
$cref fix_likelihood$$ is empty,
$code fix_like_jac_$$ is not modified.
Upon return, $code fix_like_jac_$$ contains the
$cref sparse_jac_rcv$$ structure for the
Jacobian corresponding to
$latex g_\theta ( \theta )$$ see
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.

$subhead fix_like_fun_$$
This ADFun object can be used, with $code fix_like_jac_$$,
for computing sparse Jacobians; see
$cref/f/sparse_jac_rcv/Computing Sparse Jacobians/f/$$.

$head fix_like_hes_$$
The input value of
$codei%
	CppAD::mixed::sparse_hes_rcv fix_like_hes_
%$$
must be empty.
If the return value for
$cref fix_likelihood$$ is empty,
$code fix_like_hes_$$ is not modified.
Upon return, $code fix_like_hes_$$ contains
$cref sparse_hes_rcv$$ for the
lower triangle of a Hessian corresponding to
$latex g_{\theta,\theta}) ( \theta )$$ see
$cref/g(theta)/theory/Fixed Likelihood, g(theta)/$$.
If $icode quasi_fixed$$ is true,
this is not used by $cref optimize_fixed$$, but it may be used by
$cref information_mat$$.

$subhead fix_like_fun_$$
This ADFun object can be used, with $code fix_like_hes_$$,
for computing sparse Hessians; see
$cref/f/sparse_hes_rcv/Computing Sparse Hessians/f/$$.

$end
*/


void cppad_mixed::init_fix_like(const d_vector& fixed_vec  )
{	assert( fixed_vec.size() == n_fixed_ );
	assert( ! init_fix_like_done_ );

	// ------------------------------------------------------------------------
	// fix_like_fun_
	// ------------------------------------------------------------------------
	// convert to an a1_vector
	a1_vector a1_theta(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		a1_theta[j] = fixed_vec[j];

	// start recording a1_double operations
	Independent(a1_theta);

	// compute fix_likelihood
	a1_vector a1_vec = fix_likelihood(a1_theta);
	if( a1_vec.size() == 0 )
	{	CppAD::AD<double>::abort_recording();
		init_fix_like_done_ = true;
		assert( fix_like_fun_.size_var() == 0 );
		return;
	}

	// save the recording
	fix_like_fun_.Dependent(a1_theta, a1_vec);

	// optimize the recording
# if CPPAD_MIXED_OPTIMIZE_CPPAD_FUNCTION
	fix_like_fun_.optimize();
# endif

	// ------------------------------------------------------------------------
	// fix_like_jac_
	// ------------------------------------------------------------------------
	// compute the sparsity pattern for the Jacobian
	sparse_rc pattern_in(n_fixed_, n_fixed_, n_fixed_);
	for(size_t k = 0; k < n_fixed_; k++)
		pattern_in.set(k, k, k);
	bool      transpose     = false;
	bool      dependency    = false;
	bool      internal_bool = bool_sparsity_;
	sparse_rc jac_pattern;
	fix_like_fun_.for_jac_sparsity(
		pattern_in, transpose, dependency, internal_bool, jac_pattern
	);

	// compute entire Jacobian
	fix_like_jac_.subset = sparse_rcv( jac_pattern );

	// use reverse mode for this sparse Jacobian
	// (expect fix_like_fun_.Range() <= fix_like_fun_.Domain())
	fix_like_jac_.forward = false;

	// use clear to make it clear that work is being computed
	// (should already be empty).
	fix_like_jac_.work.clear();

	// compute the work vector for reuse during sparse Jacobian calculations
	std::string coloring = "cppad";
	fix_like_fun_.sparse_jac_rev(
		fixed_vec            ,
		fix_like_jac_.subset ,
		jac_pattern          ,
		coloring             ,
		fix_like_jac_.work
	);
	// ------------------------------------------------------------------------
	// fix_like_hes_
	// ------------------------------------------------------------------------
	// no need to recalculate forward sparsity pattern.
	// (already have forward sparstity pattern stored)
	//
	// sparsity pattern for the Hessian
	size_t m = fix_like_fun_.Range();
	CppAD::vector<bool> select_range(m);
	for(size_t i = 0; i < m; i++)
		select_range[i] = false;
	select_range[0] = true;
	sparse_rc hes_pattern;
	fix_like_fun_.rev_hes_sparsity(
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
	fix_like_hes_.subset = sparse_rcv( lower_pattern );

	// compute work for reuse during sparse Hessian calculations
	// (value of weight vector does not affect work, but avoid warnings)
	d_vector weight(m);
	for(size_t i = 0; i < m; i++)
		weight[i] = 1.0;
	coloring  = "cppad.symmetric";
	fix_like_fun_.sparse_hes(
		fixed_vec            ,
		weight               ,
		fix_like_hes_.subset ,
		hes_pattern          ,
		coloring             ,
		fix_like_hes_.work
	);
	//
	init_fix_like_done_ = true;
	return;
}
