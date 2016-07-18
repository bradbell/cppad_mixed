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
$begin init_ran_objcon_hes$$
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
$icode%mixed_object%.init_ran_objcon_hes(%fixed_vec%, %random_vec%)%$$


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


$head ran_objcon_hes_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info ran_objcon_hes_
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

namespace { // BEGIN_EMPTY_NAMESPACE

// --------------------------------------------------------------------------
// ran_objcon_hes_use_set
// --------------------------------------------------------------------------
void ran_objcon_hes_use_set(
	size_t                             n_fixed      ,
	size_t                             n_random     ,
	const CppAD::vector<double>&       beta_theta_u ,
	CppAD::ADFun<double>&              fun          ,
	CppAD::mixed::sparse_hes_info&     hes_info     )
{	typedef CppAD::vector< std::set<size_t> > set_sparsity;
	hes_info.work.clear();
	assert( hes_info.row.size() == 0 );
	assert( hes_info.col.size() == 0 );
	//
	// total number of variables in H(beta, theta, u)
	size_t n_total = 2 * n_fixed + n_random;
	assert( beta_theta_u.size() == n_total );
	// -----------------------------------------------------------------------
	// compute Jacobian sparsity corresponding to partials w.r.t beta
	set_sparsity r(n_total);
	for(size_t i = 0; i < n_fixed; i++)
		r[i].insert(i);
	fun.ForSparseJac(n_fixed, r);
	// -----------------------------------------------------------------------
	// compute sparsity pattern corresponding to
	// partial w.r.t. (beta, theta, u) of partial w.r.t beta of H
	set_sparsity s(1), pattern;
	assert( s[0].empty() );
	s[0].insert(0);
	bool transpose = true;
	pattern = fun.RevSparseHes(n_fixed, s, transpose);
	// -----------------------------------------------------------------------
	// Set the row and column indices
	//
	// subsample to just the beta indices
	std::set<size_t>::iterator itr;
	for(size_t i = 0; i < n_fixed; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++ )
		{	size_t j = *itr;
			// only include lower triangular part
			if( j <= i )
			{	hes_info.row.push_back(i);
				hes_info.col.push_back(j);
			}
		}
	}
	// -----------------------------------------------------------------------
	// create a weighting vector
	CppAD::vector<double> w( fun.Range() );
	w[0] = 1.0;
	//
	// place where results go
	CppAD::vector<double> not_used( hes_info.row.size() );
	//
	// compute the work vector
	fun.SparseHessian(
		beta_theta_u,
		w,
		pattern,
		hes_info.row,
		hes_info.col,
		not_used,
		hes_info.work
	);
	return;
}

// --------------------------------------------------------------------------
// ran_objcon_hes_use_bool
// --------------------------------------------------------------------------
void ran_objcon_hes_use_bool(
	size_t                             n_fixed      ,
	size_t                             n_random     ,
	const CppAD::vector<double>&       beta_theta_u ,
	CppAD::ADFun<double>&              fun          ,
	CppAD::mixed::sparse_hes_info&     hes_info     )
{	typedef CppAD::vectorBool bool_sparsity;
	hes_info.work.clear();
	assert( hes_info.row.size() == 0 );
	assert( hes_info.col.size() == 0 );
	//
	// total number of variables in H(beta, theta, u)
	size_t n_total = 2 * n_fixed + n_random;
	assert( beta_theta_u.size() == n_total );
	// -----------------------------------------------------------------------
	// compute Jacobian sparsity corresponding to partials w.r.t beta
	bool_sparsity r(n_total * n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
			r[i * n_fixed + j] = (i == j);
	}
	fun.ForSparseJac(n_fixed, r);
	// -----------------------------------------------------------------------
	// compute sparsity pattern corresponding to
	// partial w.r.t. (beta, theta, u) of partial w.r.t beta of H
	bool_sparsity s( fun.Range() ), pattern;
	s[0] = true;
	for(size_t j = 1; j < s.size(); j++)
		s[j] = false;
	bool transpose = true;
	pattern = fun.RevSparseHes(n_fixed, s, transpose);
	// -----------------------------------------------------------------------
	// Set the row and column indices
	//
	// subsample to just the beta indices
	for(size_t i = 0; i < n_fixed; i++)
	{	// only include lower triangular part
		for(size_t j = 0; j <= i; j++)
		{	if( pattern[ i * n_fixed + j ] )
			{	hes_info.row.push_back(i);
				hes_info.col.push_back(j);
			}
		}
	}
	// -----------------------------------------------------------------------
	// create a weighting vector
	CppAD::vector<double> w( fun.Range() );
	w[0] = 1.0;
	//
	// place where results go
	CppAD::vector<double> not_used( hes_info.row.size() );
	//
	// extend sparsity pattern to all the variables
	bool_sparsity extended_pattern(n_total * n_total);
	for(size_t i = 0; i < n_total; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
			extended_pattern[i * n_total + j] =
				pattern[ i * n_fixed + j ];
		for(size_t j = n_fixed; j < n_total; j++)
			extended_pattern[i * n_total + j] = false;
	}
	//
	// compute the work vector
	fun.SparseHessian(
		beta_theta_u,
		w,
		extended_pattern,
		hes_info.row,
		hes_info.col,
		not_used,
		hes_info.work
	);
	return;
}
} // END_EMPTY_NAMESPACE

void cppad_mixed::init_ran_objcon_hes(
	bool            bool_sparsity ,
	const d_vector& fixed_vec     ,
	const d_vector& random_vec    )
{	assert( ! init_ran_objcon_hes_done_ );
	assert( init_ran_objcon_done_ );

	// total number of variables in H
	size_t n_total = 2 * n_fixed_ + n_random_;

	//	create an a1d_vector containing (theta, theta , u)
	d_vector beta_theta_u(n_total);
	pack(fixed_vec, fixed_vec, random_vec, beta_theta_u);

	if( bool_sparsity ) ran_objcon_hes_use_bool(
		n_fixed_,
		n_random_,
		beta_theta_u,
		ran_objcon_fun_,
		ran_objcon_hes_
	);
	else ran_objcon_hes_use_set(
		n_fixed_,
		n_random_,
		beta_theta_u,
		ran_objcon_fun_,
		ran_objcon_hes_
	);
	init_ran_objcon_hes_done_ = true;
}
