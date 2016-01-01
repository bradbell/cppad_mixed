// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/chol_hes_ran.hpp>

/*
$begin init_hes_ran$$
$spell
	CppAD
	init
	cppad
	hes hes
	vec
	const
	Cpp
	logdet
	Cholesky
	namespace
	Eigen
	hpp
	Simplicial
	triangular
	chol
	dismod
$$

$section Initialize Hessian of Random Likelihood w.r.t Random Effects$$

$head Syntax$$
$icode%mixed_object%.init_hes_ran(%fixed_vec%, %random_vec%)%$$

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

$head hes_ran_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info hes_ran_
%$$
does not matter.
Upon return it contains the
$cref sparse_hes_info$$
for the lower triangle of the Hessian
$latex \[
	f_{u,u} ( \theta , u )
\]$$
see $cref/f(theta, u)/
	theory/
	Random Likelihood, f(theta, u)
/$$

$subhead ran_likelihood_fun_, ran_likelihood_a1fun_$$
Either $code ran_likelihood_fun_$$ or $code ran_likelihood_a1fun_$$
can be used for the ADFun object in the
$cref/sparse Hessian Call/sparse_hes_info/Sparse Hessian Call/f/$$.

$list number$$
The matrix is symmetric and hence can be recovered from
its lower triangle.
$lnext
These indices are relative to both the fixed and random effects
with the fixed effects coming first.
$lnext
You can replace the $code a1_double$$ vectors by $code double$$ vectors,
and replace $code ran_likelihood_a1fun_$$ by $code ran_likelihood_fun_$$,
and get the results in $code double$$ instead of $code a1_double$$.
$lend

$head Random Effects Index$$
To get the indices relative to just the random effects, subtract
$code n_fixed_$$; i.e.,
$codei%
	hes_ran_.row[%k%] - n_fixed_
	hes_ran_.col[%k%] - n_fixed_
%$$
are between zero and the $code n_random_$$ and
are the row and column indices for the Hessian element
that corresponds to the $th k$$ component of
$icode a1_val_out$$ in the call to $code SparseHessian$$.

$head Lower Triangle$$
The result are only for the lower triangle of the Hessian; i.e.,
$codei%
	hes_ran_.row[%k%] >= hes_ran_.col[%k%]
%$$

$head Order$$
The results are in column major order; i.e.,
$codei%
	hes_ran_.col[%k%] <= hes_ran_.col[%k+1%]
	if( hes_ran_.col[%k%] == hes_ran_.col[%k+1%] )
		hes_ran_.row[%k%] < hes_ran_.row[%k+1%]
%$$

$head hes_ran_fun_$$
The input value of the member variables
$codei%
	CppAD::ADFun<double> hes_ran_fun_
%$$
does not matter.
Upon return its zero order forward mode computes
the lower triangle of the sparse Hessian
$latex \[
	f_{u,u} ( \theta , u )
\]$$
in the same order as the $icode a1_val_out$$ above.

$head CppAD::mixed::chol_hes_ran_$$
The static $code CppAD::mixed::chol_hes_ran_$$ is initialized by a call to
$cref/analyze_chol_hes_ran/chol_hes_ran/analyze_chol_hes_ran/$$.

$contents%example/private/hes_ran_fun_xam.cpp
%$$

$end
*/



void cppad_mixed::init_hes_ran(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( ! init_hes_ran_done_ );
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
	// compute Jacobian sparsity corresponding to parital w.r.t. random effects
	typedef CppAD::vector< std::set<size_t> > sparsity_pattern;
	sparsity_pattern r(n_total);
	for(i = n_fixed_; i < n_total; i++)
		r[i].insert(i);
	ran_likelihood_fun_.ForSparseJac(n_total, r);

	// compute sparsity pattern corresponding to paritls w.r.t. (theta, u)
	// of partial w.r.t. u of f(theta, u)
	bool transpose = true;
	sparsity_pattern s(1), pattern;
	assert( s[0].empty() );
	s[0].insert(0);
	pattern = ran_likelihood_fun_.RevSparseHes(n_total, s, transpose);


	// determine row and column indices in lower triangle of Hessian
	// and set key for column major sorting
	CppAD::vector<size_t> row, col, key;
	std::set<size_t>::iterator itr;
	for(i = n_fixed_; i < n_total; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	j = *itr;
			assert( j >= n_fixed_ );
			// only compute lower triangular part of Hessian w.r.t u only
			if( i >= j )
			{	row.push_back(i);
				col.push_back(j);
				key.push_back( i + j * n_total );
			}
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

	// sparsity pattern for range space of function
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
		ran_likelihood_fun_.ForSparseJac(n_column, r);

		// compute sparsity pattern corresponding to paritls w.r.t. (theta, u)
		// of partial w.r.t. the selected columns
		bool transpose = true;
		s[0] = true;
		h = ran_likelihood_fun_.RevSparseHes(n_column, s, transpose);

		// fill in the corresponding columns of total_sparsity
		for(i = 0; i < n_total; i++)
		{	for(j = 0; j < n_column; j++)
			{	if( i_column + j < n_total )
					pattern[i * n_total + i_column + j] = h[i * n_column + j];
			}
		}
	}
	// determine row and column indices in lower triangle of Hessian
	// and set key for column major sorting
	CppAD::vector<size_t> row, col, key;
	for(i = n_fixed_; i < n_total; i++)
	{	for(j = n_fixed_; j < n_total; j++)
		{	if( pattern[i * n_total + j] )
			{	// only compute lower triangular of Hessian w.r.t u only
				if( i >= j )
				{	row.push_back(i);
					col.push_back(j);
					key.push_back( i + j * n_total );
				}
			}
		}
	}
# endif
	// -----------------------------------------------------------------------

	// set hes_ran_.row and hes_ran_.col in colum major order
	size_t K = row.size();
	CppAD::vector<size_t> ind(K);
	CppAD::index_sort(key, ind);
	hes_ran_.row.resize(K);
	hes_ran_.col.resize(K);
	for(size_t k = 0; k < row.size(); k++)
	{	hes_ran_.row[k] = row[ ind[k] ];
		hes_ran_.col[k] = col[ ind[k] ];
	}

	// create a weighting vector
	d_vector w(1);
	w[0] = 1.0;

	// place where results go (not usd here)
	d_vector val_out(K);

	// compute the work vector
	ran_likelihood_fun_.SparseHessian(
		both,
		w,
		pattern,
		hes_ran_.row,
		hes_ran_.col,
		val_out,
		hes_ran_.work
	);

	// now tape the same computation and store in hes_ran_fun_
	CppAD::vector< std::set<size_t> > not_used(0);
	a1d_vector a1_both(n_total), a1_w(1), a1_val_out(K);
	for(size_t i = 0; i < n_total; i++)
		a1_both[i] = both[i];
	a1_w[0] = w[0];
	CppAD::Independent(a1_both);
	ran_likelihood_a1fun_.SparseHessian(
		a1_both,
		a1_w,
		not_used,
		hes_ran_.row,
		hes_ran_.col,
		a1_val_out,
		hes_ran_.work
	);
	hes_ran_fun_.Dependent(a1_both, a1_val_out);
	//
	init_hes_ran_done_ = true;

	// now analyze the lower triangular Cholesky factorization
	CppAD::mixed::analyze_chol_hes_ran(n_fixed_, n_random_, hes_ran_.row, hes_ran_.col);
}


