// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad_mixed/cppad_mixed.hpp>

/*
$begin init_hes_ranobj$$
$spell
	init
	ranobj
	cppad
	obj
	hes
	vec
	const
	Cpp
$$

$section Initialize Hessian of Approximate Random Objective$$

$head Syntax$$
$codei%init_hes_ranobj(%fixed_vec%, %random_vec%)%$$

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

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the initial value for the
$cref/random effects/cppad_mixed/Random Effects, u/$$ optimization.


$head hes_ranobj_$$
The input value of the member variable
$codei%
	sparse_hes_info hes_ranobj_
%$$
does not matter.
Upon return it contains the
$cref sparse_hes_info$$
for the Hessian
$latex \[
	r^{(2)} ( \theta )
	=
	H_{\beta \beta}^{(2)} [ \beta, \theta , \hat{u} ( \theta) ]
\] $$
see
$cref/H(beta, theta, u)
	/cppad_mixed_theory
	/Hessian of Random Objective
	/Approximate Random Objective, H(beta, theta, u)
/$$.
Note that the matrix is symmetric and hence can be recovered from
its lower triangle.

$subhead ranobj_fun_$$
This ADFun object should be used in the
$cref/sparse Hessian Call/sparse_hes_info/Sparse Hessian Call/f/$$.

$end
*/

namespace cppad_mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void cppad_mixed::init_hes_ranobj(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_hes_ranobj_done_ );
	assert( init_ranobj_done_ );
	size_t i, j;

	// total number of variables in H
	size_t n_total = 2 * n_fixed_ + n_random_;

	//	create an a1d_vector containing (theta, theta , u)
	d_vector beta_theta_u(n_total);
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);


	// compute Jacobian sparsity corresponding to partial w.r.t beta
	// of H(beta, beta, u)
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	sparsity_pattern r(n_total);
	for(i = 0; i < n_fixed_; i++)
		r[i].insert(i);
	ranobj_fun_.ForSparseJac(n_fixed_, r);

	// compute sparsity pattern corresponding to partial w.r.t (beta, theta, u)
	// of parital w.r.t beta of H(beta, theta, u)
	sparsity_pattern s(1);
	assert( s[0].empty() );
	s[0].insert(0);
	bool transpose = true;
	sparsity_pattern pattern =
		ranobj_fun_.RevSparseHes(n_fixed_, s, transpose);

	// determine row and column indices in lower triangle of Hessian
	hes_ranobj_.row.clear();
	hes_ranobj_.col.clear();
	std::set<size_t>::iterator itr;
	for(i = 0; i < n_fixed_; i++)
	{	for(
			itr  = pattern[i].begin();
			itr != pattern[i].end();
			itr++
		)
		{	j = *itr;
			// only compute lower triangular part
			if( i >= j )
			{	hes_ranobj_.row.push_back(i);
				hes_ranobj_.col.push_back(j);
			}
		}
	}

	// create a weighting vector
	d_vector w(1);
	w[0] = 1.0;

	// place where results go (not used here)
	d_vector val_out( hes_ranobj_.row.size() );

	// compute the work vector
	ranobj_fun_.SparseHessian(
		beta_theta_u,
		w,
		pattern,
		hes_ranobj_.row,
		hes_ranobj_.col,
		val_out,
		hes_ranobj_.work
	);
	//
	init_hes_ranobj_done_ = true;
}

} // END_CPPAD_MIXED_NAMESPACE
