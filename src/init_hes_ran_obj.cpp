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
$begin init_hes_ran_obj$$
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
$$

$section Initialize Hessian of Approximate Random Objective$$

$head Syntax$$
$icode%mixed_object%.init_hes_ran_obj(%fixed_vec%, %random_vec%)%$$


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

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the initial value for the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ optimization.


$head hes_ran_obj_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info hes_ran_obj_
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
	/Approximate Random Objective, H(beta, theta, u)
/$$.
Note that the matrix is symmetric and hence can be recovered from
its lower triangle.

$subhead ran_objcon_fun_$$
This ADFun object should be used in the
$cref/sparse Hessian Call/sparse_hes_info/Sparse Hessian Call/f/$$.

$end
*/


void cppad_mixed::init_hes_ran_obj(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_hes_ran_obj_done_ );
	assert( init_ran_objcon_done_ );
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
	ran_objcon_fun_.ForSparseJac(n_fixed_, r);

	// compute sparsity pattern corresponding to partial w.r.t (beta, theta, u)
	// of parital w.r.t beta of H(beta, theta, u)
	sparsity_pattern s(1);
	assert( s[0].empty() );
	s[0].insert(0);
	bool transpose = true;
	sparsity_pattern pattern =
		ran_objcon_fun_.RevSparseHes(n_fixed_, s, transpose);

	// determine row and column indices in lower triangle of Hessian
	hes_ran_obj_.row.clear();
	hes_ran_obj_.col.clear();
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
			{	hes_ran_obj_.row.push_back(i);
				hes_ran_obj_.col.push_back(j);
			}
		}
	}

	// create a weighting vector
	d_vector w(1);
	w[0] = 1.0;

	// place where results go (not used here)
	d_vector val_out( hes_ran_obj_.row.size() );

	// compute the work vector
	ran_objcon_fun_.SparseHessian(
		beta_theta_u,
		w,
		pattern,
		hes_ran_obj_.row,
		hes_ran_obj_.col,
		val_out,
		hes_ran_obj_.work
	);
	//
	init_hes_ran_obj_done_ = true;
}

