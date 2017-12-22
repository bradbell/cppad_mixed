// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin init_laplace_obj_hes$$
$spell
	objcon
	CppAD
	init
	ran_obj
	cppad
	obj
	hes
	vec
	const
	Cpp
	bool
$$

$section Initialize Hessian of Approximate Laplace Objective$$

$head Syntax$$
$icode%mixed_object%.init_laplace_obj_hes(
	%fixed_vec%, %random_vec%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Assumptions$$
The member variable
$cref/init_laplace_obj_done_/init_ran_objcon/init_laplace_obj_done_/$$
is true.

$head init_laplace_obj_hes_done_$$
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

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the initial value for the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ optimization.


$head laplace_obj_hes_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info laplace_obj_hes_
%$$
does not matter.
Upon return it contains the
$cref sparse_hes_info$$
for the lower triangle of the Hessian
$latex \[
	r_{\theta,\theta} ( \theta )
	=
	H_{\beta,\beta} [ \beta, \theta , \hat{u} ( \theta) ]
\] $$
see
$cref/H(beta, theta, u)
	/theory
	/Approximate Laplace Objective, H(beta, theta, u)
/$$.
Note that the matrix is symmetric and hence can be recovered from
its lower triangle.

$subhead laplace_obj_fun_$$
This ADFun object should be used in the
$cref/sparse Hessian Call/sparse_hes_info/Sparse Hessian Call/f/$$.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_laplace_obj_hes(
	const d_vector& fixed_vec     ,
	const d_vector& random_vec    )
{	assert( ! init_laplace_obj_hes_done_ );
	assert( init_laplace_obj_done_ );

	// total number of variables in H
	size_t n_total = 2 * n_fixed_ + n_random_;

	//	create an a1_vector containing (theta, theta , u)
	d_vector beta_theta_u(n_total);
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);

	// compute Jacobian sparsity for partials w.r.t. beta
	sparse_rc pattern_in(n_total, n_fixed_, n_fixed_);
	for(size_t k = 0; k < n_fixed_; k++)
		pattern_in.set(k, k, k);
	bool      transpose     = false;
	bool      dependency    = false;
	bool      internal_bool = bool_sparsity_;
	sparse_rc jac_pattern;
	laplace_obj_fun_.for_jac_sparsity(
		pattern_in, transpose, dependency, internal_bool, jac_pattern
	);

	// compute sparsity pattern corresponding to partial w.r.t (beta, theta, u)
	// of partial w.r.t beta of H(beta, theta, u)
	size_t m = laplace_obj_fun_.Range();
	CppAD::vector<bool> select_range(m);
	for(size_t i = 0; i < m; i++)
		select_range[i] = false;
	select_range[0] = true;
	sparse_rc hes_pattern;
	laplace_obj_fun_.rev_hes_sparsity(
		select_range, transpose, internal_bool, hes_pattern
	);
	assert( hes_pattern.nr() == n_total );
	assert( hes_pattern.nc() == n_fixed_ );
	//
	// sparsity pattern corresponding to lower traingle of
	// partial w.r.t. beta of partial w.r.t. beta of H(beta, theta, u)
	size_t nnz = 0;
	for(size_t k = 0; k < hes_pattern.nnz(); k++)
	{	size_t r = hes_pattern.row()[k];
		size_t c = hes_pattern.col()[k];
		if( n_fixed_ > r && r >= c )
			++nnz;
	}
	sparse_rc beta_beta_pattern(n_total, n_total, nnz);
	sparse_rc extend_pattern(n_total, n_total, hes_pattern.nnz());
	size_t ell = 0;
	for(size_t k = 0; k < hes_pattern.nnz(); k++)
	{	size_t r = hes_pattern.row()[k];
		size_t c = hes_pattern.col()[k];
		extend_pattern.set(k, r, c);
		if( n_fixed_ > r && r >= c )
			beta_beta_pattern.set(ell++, r, c);
	}
	assert( nnz == ell );
	//
	// only compute the lower triangle of H_beta_beta ( beta, theta, u)
	laplace_obj_hes_.subset = sparse_rcv( beta_beta_pattern );
	//
	// compute work for reuse during sparse Hessian calculations
	// (value of weight vecgttor does not affectg work, but avoid wanrings)
	d_vector weight(m);
	for(size_t i = 0; i < m; i++)
		weight[i] = 1.0;
	std::string coloring  = "cppad.symmetric";
	laplace_obj_fun_.sparse_hes(
		beta_theta_u           ,
		weight                 ,
		laplace_obj_hes_.subset ,
		extend_pattern         ,
		coloring               ,
		laplace_obj_hes_.work
	);
	//
	init_laplace_obj_hes_done_ = true;
}
