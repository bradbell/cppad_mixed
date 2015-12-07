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
# include <dismod_at/configure.hpp>

/*
$begin init_hes_cross$$
$spell
	init
	cppad
	hes hes
	vec
	const
	Cpp
	logdet
$$

$section Cross Terms of Sparse Hessian w.r.t Fixed and Random Effects$$

$head Syntax$$
$codei%init_hes_cross(%fixed_vec%, %random_vec%)%$$

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
It specifies the value of the
$cref/random effects/cppad_mixed/Random Effects, u/$$
vector $latex u$$ at which the initialization is done.

$head hes_cross_$$
The input value of the member variable
$codei%
	sparse_hes_info hes_cross_
%$$
does not matter.
Upon return it contains the
$cref sparse_hes_info$$
for the Hessian
$latex \[
	f_{u \theta}^{(2)} ( \theta , u )
\]$$
see $cref/f(theta, u)/
	cppad_mixed_theory/
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
	example/devel/cppad_mixed/private/hes_cross_xam.cpp
%$$
$head Example$$
The file $cref hes_cross_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


namespace dismod_at { // BEGIN_DISMOD_AT_NAMESPACE

void cppad_mixed::init_hes_cross(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{
	assert( ! init_hes_cross_done_ );
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	size_t i, j;

	// total number of variables
	size_t n_total = n_fixed_ + n_random_;

	// create a d_vector containing (theta, u)
	d_vector both(n_total);
	pack(fixed_vec, random_vec, both);

# if CPPAD_MIXED_SET_SPARSITY
	// ----------------------------------------------------------------------
	// compute Jacobian sparsity corresponding to parital w.r.t. fixed effects
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	sparsity_pattern r(n_total);
	for(i = 0; i < n_fixed_; i++)
		r[i].insert(i);
	ran_like_fun_.ForSparseJac(n_total, r);

	// compute sparsity pattern corresponding to paritls w.r.t. (theta, u)
	// of partial w.r.t. theta of f(theta, u)
	bool transpose = true;
	sparsity_pattern s(1), pattern;
	assert( s[0].empty() );
	s[0].insert(0);
	pattern = ran_like_fun_.RevSparseHes(n_total, s, transpose);


	// User row index for random effect and column index for fixed effect
	CppAD::vector<size_t> row, col, key;
	std::set<size_t>::iterator itr;
	for(i = n_fixed_; i < n_total; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	j = *itr;
			assert( j < n_fixed_ );
			row.push_back(i);
			col.push_back(j);
			key.push_back( i + j * n_total );
		}
	}
# else
	// -----------------------------------------------------------------------
	// compute Jacobian sparsity corresponding to parital w.r.t. random effects
	typedef CppAD::vectorBool sparsity_pattern;

	// sparsity pattern for complete Hessian
	sparsity_pattern pattern(n_total * n_total);

	// number of bits that are packed into one unit in vectorBool
	// (note yet part of CppAD API).
	size_t n_column = std::numeric_limits<size_t>::digits;

	// sparsity patterns for current columns
	sparsity_pattern r(n_total * n_column), h(n_total * n_column);

	// sparsity patter for range space of function
	sparsity_pattern s(1);

	// loop that computes the sparsity pattern n_column columns at a time
	size_t n_loop = (n_total - 1) / n_column + 1;
	for(size_t i_loop = 0; i_loop < n_loop; i_loop++)
	{	// starting column index for this iteration
		size_t i_column = i_loop * n_column;

		// pattern that picks out the appropriate columns
		for(i = 0; i < n_total; i++)
		{	for(j = 0; j < n_column; j++)
				r[i * n_column + j] = (i == i_column + j);
		}
		ran_like_fun_.ForSparseJac(n_column, r);

		// compute sparsity pattern corresponding to paritls w.r.t. (theta, u)
		// of partial w.r.t. the selected columns
		bool transpose = true;
		s[0] = true;
		h = ran_like_fun_.RevSparseHes(n_column, s, transpose);

		// fill in the corresponding columns of total_sparsity
		for(i = 0; i < n_total; i++)
		{	for(j = 0; j < n_column; j++)
			{	if( i_column + j < n_total )
					pattern[i * n_total + i_column + j] = h[i * n_column + j];
			}
		}
	}

	// User row index for random effect and column index for fixed effect
	CppAD::vector<size_t> row, col, key;
	for(i = n_fixed_; i < n_total; i++)
	{	for(j = 0; j < n_fixed_; j++)
		{	if( pattern[i * n_total + j] )
			{	row.push_back(i);
				col.push_back(j);
				key.push_back( i + j * n_total );
			}
		}
	}
# endif
	// ----------------------------------------------------------------------
	// set hes_cross_.row and hes_cross_.col in colum major order
	size_t K = row.size();
	CppAD::vector<size_t> ind(K);
	CppAD::index_sort(key, ind);
	hes_cross_.row.resize(K);
	hes_cross_.col.resize(K);
	for(size_t k = 0; k < row.size(); k++)
	{	hes_cross_.row[k] = row[ ind[k] ];
		hes_cross_.col[k] = col[ ind[k] ];
	}

	// create a weighting vector
	d_vector w(1);
	w[0] = 1.0;

	// place where results go (not usd here)
	d_vector val_out(K);

	// compute the work vector
	ran_like_fun_.SparseHessian(
		both,
		w,
		pattern,
		hes_cross_.row,
		hes_cross_.col,
		val_out,
		hes_cross_.work
	);
	//
	init_hes_cross_done_ = true;
}


} // END_DISMOD_AT_NAMESPACE
