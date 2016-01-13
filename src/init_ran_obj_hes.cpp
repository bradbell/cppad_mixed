// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

/*
$begin init_ran_obj_hes$$
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
$icode%mixed_object%.init_ran_obj_hes(%fixed_vec%, %random_vec%)%$$


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


$head ran_obj_hes_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info ran_obj_hes_
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
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/configure.hpp>


void cppad_mixed::init_ran_obj_hes(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_ran_obj_hes_done_ );
	assert( init_ran_objcon_done_ );
	assert( init_ran_con_done_ );
	size_t i, j;

	// total number of variables in H
	size_t n_total = 2 * n_fixed_ + n_random_;

	//	create an a1d_vector containing (theta, theta , u)
	d_vector beta_theta_u(n_total);
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);

	// 2DO: use CPPAD_MIXED_SET_SPARSITY to choose between set
	// and vectorBool sparsity

	// compute Jacobian sparsity corresponding to partial w.r.t beta
	// of [ H(beta, theta, u) , B(beta, theta, u) ]
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	sparsity_pattern r(n_total);
	for(i = 0; i < n_fixed_; i++)
		r[i].insert(i);
	ran_objcon_fun_.ForSparseJac(n_fixed_, r);

	// compute sparsity pattern corresponding to partial w.r.t (beta, theta, u)
	// of parital w.r.t beta of H(beta, theta, u) which is the first component
	// of ran_obj_fun_.
	sparsity_pattern s(1);
	assert( s[0].empty() );
	s[0].insert(0);
	bool transpose = true;
	sparsity_pattern pattern =
		ran_objcon_fun_.RevSparseHes(n_fixed_, s, transpose);

	// determine row and column indices in lower triangle of Hessian
	ran_obj_hes_.row.clear();
	ran_obj_hes_.col.clear();
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
			{	ran_obj_hes_.row.push_back(i);
				ran_obj_hes_.col.push_back(j);
			}
		}
	}

	// create a weighting vector
	d_vector w(1 + n_ran_con_);
	w[0] = 1.0;

	// place where results go (not used here)
	d_vector val_out( ran_obj_hes_.row.size() );

	// compute the work vector
	ran_objcon_fun_.SparseHessian(
		beta_theta_u,
		w,
		pattern,
		ran_obj_hes_.row,
		ran_obj_hes_.col,
		val_out,
		ran_obj_hes_.work
	);
	//
	init_ran_obj_hes_done_ = true;
}
