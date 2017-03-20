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
$begin init_hes_cross$$
$spell
	CppAD
	init
	cppad
	hes hes
	vec
	const
	Cpp
	logdet
	bool
$$

$section Cross Terms of Sparse Hessian w.r.t Fixed and Random Effects$$

$head Syntax$$
$icode%mixed_object%.init_hes_cross(
	%fixed_vec%, %random_vec%
)%$$

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
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which the initialization is done.

$head hes_cross_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info hes_cross_
%$$
does not matter.
Upon return it contains the
$cref sparse_hes_info$$
for the Hessian
$latex \[
	f_{u,\theta} ( \theta , u )
\]$$
see $cref/f(theta, u)/
	theory/
	Random Likelihood, f(theta, u)
/$$

$subhead ran_like_fun_, ran_like_a1fun_$$
Either $code ran_like_fun_$$ or $code ran_like_a1fun_$$
can be used for the ADFun object in the
$cref/sparse Hessian Call/sparse_hes_info/Sparse Hessian Call/f/$$.


$head Order$$
The results are in column major order; i.e.,
$codei%
	hes_cross_.col[%k%] <= hes_cross_.col[%k+1%]
	if( hes_cross_.col[%k%] == hes_cross_.col[%k+1%] )
		hes_cross_.row[%k%] < hes_cross_.row[%k+1%]
%$$

$children%
	example/private/hes_cross.cpp
%$$
$head Example$$
The file $cref hes_cross.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/

namespace { // BEGIN_EMPTY_NAMESPACE
// --------------------------------------------------------------------------
// cross_hes_use_set
// --------------------------------------------------------------------------
void cross_hes_use_set(
	size_t                             n_fixed  ,
	size_t                             n_random ,
	const CppAD::vector<double>&       both     ,
	CppAD::ADFun<double>&              fun      ,
	CppAD::mixed::sparse_hes_info&     hes_info )
{	hes_info.work.clear();
	assert( hes_info.row.size() == 0 );
	assert( hes_info.col.size() == 0 );
	assert( fun.Range() == 1 );
	assert( both.size() == n_fixed + n_random );
	// ----------------------------------------------------------------------
	typedef CppAD::vector< std::set<size_t> > set_sparsity;
	size_t n_both = n_fixed + n_random;
	// ----------------------------------------------------------------------
	// compute Jacobian sparsity corresponding to parital w.r.t. fixed effects
	set_sparsity r(n_both);
	for(size_t i = 0; i < n_fixed; i++)
		r[i].insert(i);
	fun.ForSparseJac(n_fixed, r);
	// ----------------------------------------------------------------------
	// compute sparsity pattern corresponding to
	// partial w.r.t. (theta, u) of partial w.r.t. theta
	bool transpose = true;
	set_sparsity s(1), pattern;
	assert( s[0].empty() );
	s[0].insert(0);
	pattern = fun.RevSparseHes(n_fixed, s, transpose);
	// ----------------------------------------------------------------------
	// Use row index for random effect and column index for fixed effect
	CppAD::vector<size_t> row, col, key;
	std::set<size_t>::iterator itr;
	//
	// subsample to just the random effects
	for(size_t i = n_fixed; i < n_both; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			//
			// partial w.r.t u[i] of partial w.r.t theta[j]
			assert( j < n_fixed );
			row.push_back(i);
			col.push_back(j);
			key.push_back( i + j * n_both );
		}
	}
	// ----------------------------------------------------------------------
	// set hes_info in colum major order
	size_t K = row.size();
	CppAD::vector<size_t> ind(K);
	CppAD::index_sort(key, ind);
	hes_info.row.resize(K);
	hes_info.col.resize(K);
	for(size_t k = 0; k < row.size(); k++)
	{	hes_info.row[k] = row[ ind[k] ];
		hes_info.col[k] = col[ ind[k] ];
	}
	// ----------------------------------------------------------------------
	// Compute work for reuse during future calls to SparseHessian
	CppAD::vector<double> w(1);
	w[0] = 1.0;
	CppAD::vector<double> not_used(K);
	// compute the work vector
	fun.SparseHessian(
		both,
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
// cross_hes_use_bool
// --------------------------------------------------------------------------
void cross_hes_use_bool(
	size_t                             n_fixed  ,
	size_t                             n_random ,
	const CppAD::vector<double>&       both     ,
	CppAD::ADFun<double>&              fun      ,
	CppAD::mixed::sparse_hes_info&     hes_info )
{	hes_info.work.clear();
	assert( hes_info.row.size() == 0 );
	assert( hes_info.col.size() == 0 );
	assert( fun.Range() == 1 );
	assert( both.size() == n_fixed + n_random );
	// ----------------------------------------------------------------------
	typedef CppAD::vectorBool bool_sparsity;
	size_t n_both = n_fixed + n_random;
	// ----------------------------------------------------------------------
	// compute Jacobian sparsity corresponding to parital w.r.t. fixed effects
	bool_sparsity r(n_both * n_fixed);
	for(size_t i = 0; i < n_fixed; i++)
	{	for(size_t j = 0; j < n_random; j++)
			r[i * n_fixed + j] = (i == j);
	}
	fun.ForSparseJac(n_fixed, r);
	// ----------------------------------------------------------------------
	// compute sparsity pattern corresponding to
	// partial w.r.t. (theta, u) of partial w.r.t. theta
	bool transpose = true;
	bool_sparsity s(1), pattern;
	s[0] = true;
	pattern = fun.RevSparseHes(n_fixed, s, transpose);
	// ----------------------------------------------------------------------
	// Use row index for random effect and column index for fixed effect
	CppAD::vector<size_t> row, col, key;
	std::set<size_t>::iterator itr;
	//
	// subsample to just the random effects
	for(size_t i = n_fixed; i < n_both; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
		{	if( pattern[ i * n_fixed + j ] )
			{	// partial w.r.t u[i] of partial w.r.t theta[j]
				row.push_back(i);
				col.push_back(j);
				key.push_back( i + j * n_both );
			}
		}
	}
	// ----------------------------------------------------------------------
	// set hes_info in colum major order
	size_t K = row.size();
	CppAD::vector<size_t> ind(K);
	CppAD::index_sort(key, ind);
	hes_info.row.resize(K);
	hes_info.col.resize(K);
	for(size_t k = 0; k < row.size(); k++)
	{	hes_info.row[k] = row[ ind[k] ];
		hes_info.col[k] = col[ ind[k] ];
	}
	// ----------------------------------------------------------------------
	// Compute work for reuse during future calls to SparseHessian
	CppAD::vector<double> w(1);
	w[0] = 1.0;
	CppAD::vector<double> not_used(K);
	//
	// extend sparsity pattern to all the variables
	bool_sparsity extended_pattern(n_both * n_both);
	for(size_t i = 0; i < n_both; i++)
	{	for(size_t j = 0; j < n_fixed; j++)
			extended_pattern[ i * n_both + j ] = pattern[ i * n_fixed + j ];
		for(size_t j = n_fixed; j < n_both; j++)
			extended_pattern[ i * n_both + j ] = false;
	}
	fun.SparseHessian(
		both,
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



void cppad_mixed::init_hes_cross(
	const d_vector& fixed_vec         ,
	const d_vector& random_vec        )
{
	assert( ! init_hes_cross_done_ );
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );

	// total number of variables
	size_t n_both = n_fixed_ + n_random_;

	// create a d_vector containing (theta, u)
	d_vector both(n_both);
	pack(fixed_vec, random_vec, both);

	if( bool_sparsity_ ) cross_hes_use_bool(
		n_fixed_,
		n_random_,
		both,
		ran_like_fun_,
		hes_cross_
	);
	else cross_hes_use_set(
		n_fixed_,
		n_random_,
		both,
		ran_like_fun_,
		hes_cross_
	);
	init_hes_cross_done_ = true;
}


