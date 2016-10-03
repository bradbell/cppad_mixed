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
$begin init_ran_hes$$
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
	hpp
	Simplicial
	triangular
	chol
	dismod
	bool
$$

$section Initialize Hessian of Random Likelihood w.r.t Random Effects$$

$head Syntax$$
$icode%mixed_object%.init_ran_hes(
	%bool_sparsity%, %fixed_vec%, %random_vec%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head bool_sparsity$$
This argument has prototype
$codei%
	bool %bool_sparsity%
%$$
If it is true, boolean sparsity patterns are used for this computation,
otherwise set sparsity patterns are used.

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

$head ran_hes_$$
The input value of the member variable
$codei%
	CppAD::mixed::sparse_hes_info ran_hes_
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

$subhead ran_like_fun_, ran_like_a1fun_$$
Either $code ran_like_fun_$$ or $code ran_like_a1fun_$$
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
and replace $code ran_like_a1fun_$$ by $code ran_like_fun_$$,
and get the results in $code double$$ instead of $code a1_double$$.
$lend

$head Random Effects Index$$
To get the indices relative to just the random effects, subtract
$code n_fixed_$$; i.e.,
$codei%
	ran_hes_.row[%k%] - n_fixed_
	ran_hes_.col[%k%] - n_fixed_
%$$
are between zero and the $code n_random_$$ and
are the row and column indices for the Hessian element
that corresponds to the $th k$$ component of
$icode a1_val_out$$ in the call to $code SparseHessian$$.

$head Lower Triangle$$
The result are only for the lower triangle of the Hessian; i.e.,
$codei%
	ran_hes_.row[%k%] >= ran_hes_.col[%k%]
%$$

$head Order$$
The results are in column major order; i.e.,
$codei%
	ran_hes_.col[%k%] <= ran_hes_.col[%k+1%]
	if( ran_hes_.col[%k%] == ran_hes_.col[%k+1%] )
		ran_hes_.row[%k%] < ran_hes_.row[%k+1%]
%$$

$head ran_hes_fun_$$
The input value of the member variables
$codei%
	CppAD::ADFun<double> ran_hes_fun_
%$$
does not matter.
Upon return its zero order forward mode computes
the lower triangle of the sparse Hessian
$latex \[
	f_{u,u} ( \theta , u )
\]$$
in the same order as the $icode a1_val_out$$ above.

$contents%example/private/ran_hes_fun.cpp
%$$

$end
*/


namespace { // BEGIN_EMPTY_NAMESPACE
// --------------------------------------------------------------------------
// ran_hes_use_set
// --------------------------------------------------------------------------
void ran_hes_use_set(
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
	// compute Jacobian sparsity corresponding to parital w.r.t. random effects
	set_sparsity r(n_both);
	for(size_t i = n_fixed; i < n_both; i++)
		r[i].insert(i);
	fun.ForSparseJac(n_both, r);
	// ----------------------------------------------------------------------
	// compute sparsity pattern corresponding to
	// partial w.r.t. (theta, u) of partial w.r.t. u
	bool transpose = true;
	set_sparsity s(1), pattern;
	assert( s[0].empty() );
	s[0].insert(0);
	pattern = fun.RevSparseHes(n_both, s, transpose);
	// ----------------------------------------------------------------------
	// Set the row and column indices
	//
	// subsample to just the random effects
	CppAD::vector<size_t> row, col, key;
	std::set<size_t>::iterator itr;
	for(size_t i = n_fixed; i < n_both; i++)
	{	for(itr = pattern[i].begin(); itr != pattern[i].end(); itr++)
		{	size_t j = *itr;
			//
			// partial w.r.t u[i] of partial w.r.t u[j]
			assert( n_fixed <= j );
			// only store lower triangle of symmetric Hessian
			if( i >= j )
			{	row.push_back(i);
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
// ran_hes_use_bool
// --------------------------------------------------------------------------
void ran_hes_use_bool(
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
	// compute Jacobian sparsity corresponding to parital w.r.t. random effects
	bool_sparsity r(n_both * n_random);
	for(size_t i = 0; i < n_both; i++)
	{	for(size_t j = 0; j < n_random; j++)
			r[i * n_random + j] = (i >= n_fixed) & ((i - n_fixed) == j);
	}
	fun.ForSparseJac(n_random, r);
	// ----------------------------------------------------------------------
	// compute sparsity pattern corresponding to
	// partial w.r.t. (theta, u) of partial w.r.t. theta
	bool transpose = true;
	bool_sparsity s(1), pattern;
	s[0] = true;
	pattern = fun.RevSparseHes(n_random, s, transpose);
	// ----------------------------------------------------------------------
	// Set the row and column indices
	//
	// subsample to just the random effects
	CppAD::vector<size_t> row, col, key;
	for(size_t i = n_fixed; i < n_both; i++)
	{	for(size_t j = 0; j < n_random; j++)
		{	if( pattern[ i * n_random + j ] & ( i >= j + n_fixed ) )
			{	// partial w.r.t u[i] of partial w.r.t theta[j]
				row.push_back(i);
				col.push_back(j + n_fixed);
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
			extended_pattern[ i * n_both + j ] = false;
		for(size_t j = 0; j < n_random; j++)
			extended_pattern[i * n_both + j + n_fixed ] =
				pattern[ i * n_random + j ];
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


void cppad_mixed::init_ran_hes(
	bool            bool_sparsity ,
	const d_vector& fixed_vec     ,
	const d_vector& random_vec    )
{	assert( ! init_ran_hes_done_ );
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );

	// total number of variables
	size_t n_total = n_fixed_ + n_random_;

	// create a d_vector containing (theta, u)
	d_vector both(n_total);
	pack(fixed_vec, random_vec, both);

	if( bool_sparsity ) ran_hes_use_bool(
		n_fixed_,
		n_random_,
		both,
		ran_like_fun_,
		ran_hes_
	);
	else ran_hes_use_set(
		n_fixed_,
		n_random_,
		both,
		ran_like_fun_,
		ran_hes_
	);

	// Declare the independent and dependent variables for taping calculation
	// of Hessian of the random likelihood w.r.t. the random effects
	a1d_vector a1_both(n_total);
	for(size_t i = 0; i < n_total; i++)
		a1_both[i] = both[i];
	CppAD::Independent(a1_both);

	// unpack the independent variables into a1_fixed and a1_random
	a1d_vector a1_fixed(n_fixed_), a1_random(n_random_);
	unpack(a1_fixed, a1_random, a1_both);

	// Evaluate ran_likelihood_hes. Its row indices are relative
	// to just the random effects, not both fixed and random effects.
	size_t K = ran_hes_.row.size();
	CppAD::vector<size_t> row(K), col(K);
	for(size_t k = 0; k < K; k++)
	{	assert( ran_hes_.row[k] >= n_fixed_ );
		assert( ran_hes_.col[k] >= n_fixed_ );
		row[k] = ran_hes_.row[k] - n_fixed_;
		col[k] = ran_hes_.col[k] - n_fixed_;
	}
	a1d_vector a1_val_out = ran_likelihood_hes(
		a1_fixed, a1_random, row, col
	);
	if( a1_val_out.size() == 0 )
	{	// The user has not defined ran_likelihood_hes, so use AD to calcuate
		// the Hessian of the randome likelihood w.r.t the random effects.
		CppAD::vector< std::set<size_t> > not_used(0);
		a1d_vector a1_w(1);
		a1_w[0] = 1.0;
		a1_val_out.resize(K);
		ran_like_a1fun_.SparseHessian(
			a1_both,
			a1_w,
			not_used,
			ran_hes_.row,
			ran_hes_.col,
			a1_val_out,
			ran_hes_.work
		);
	}
	ran_hes_fun_.Dependent(a1_both, a1_val_out);
	//
	init_ran_hes_done_ = true;
}

/*
$begin init_ran_hes_check$$
$spell
	cppad
	const
	CppAD
	std
	vec
	hes
	init
$$

$section Check User Defined ran_likelihood_hes$$

$head Syntax$$
$icode%mixed_object%.init_ran_hes_check(%fixed_vec%, %random_vec%)%$$

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
vector $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.

$head ran_like_fun_$$
The member variable
$codei%
	CppAD::ADFun<double> ran_like_fun_
%$$
must contain a recording of the random likelihood function; see
$cref/ran_like_fun_/init_ran_like/ran_like_fun_/$$.

$head ran_hes_$$
The member variable
$codei%
	CppAD::mixed::sparse_hes_info ran_hes_
%$$
must contain the information for evaluating the random effects
Hessian using $code ran_like_fun_$$; see
$cref/ran_hes_/init_ran_hes/ran_hes_/$$.

$head ran_likelihood_hes$$
If the return value from $cref ran_likelihood_hes$$ is non-empty,
it is checked for correctness.
To be specific, the return vector is checked using
$cref/ran_like_fun_/Private/ran_like_fun_/$$.
If the results do no agree,
$cref/fatal_error/Public/User Defined Functions/fatal_error/$$ is called with
an appropriate error message.

$end
-----------------------------------------------------------------------------
*/
void cppad_mixed::init_ran_hes_check(
	const d_vector&        fixed_vec       ,
	const d_vector&        random_vec      )
{
	// a1_double versions of fixed_vec and random_vec
	a1d_vector a1_fixed(n_fixed_), a1_random(n_random_);
	for(size_t i = 0; i < n_fixed_; i++)
		a1_fixed[i] = fixed_vec[i];
	for(size_t i = 0; i < n_random_; i++)
		a1_random[i] = random_vec[i];

	// row indices in ran_likelihood_hes are relative
	// to just the random effects, not both fixed and random effects.
	size_t K = ran_hes_.row.size();
	CppAD::vector<size_t> row(K), col(K);
	for(size_t k = 0; k < K; k++)
	{	assert( ran_hes_.row[k] >= n_fixed_ );
		assert( ran_hes_.col[k] >= n_fixed_ );
		row[k] = ran_hes_.row[k] - n_fixed_;
		col[k] = ran_hes_.col[k] - n_fixed_;
	}
	a1d_vector a1_val_out = ran_likelihood_hes(
		a1_fixed, a1_random, row, col
	);
	if( a1_val_out.size() == 0 )
		return;
	//
	// same order as chose by cppad/mixed/pack.hpp
	d_vector both( n_fixed_ + n_random_ );
	pack(fixed_vec, random_vec, both);
	//
	// weighting vector
	assert( ran_like_fun_.Range() == 1 );
	d_vector w(1);
	w[0] = 1.0;
	//
	// sparsity pattern is not used
	CppAD::vector< std::set<size_t> > not_used(0);
	//
	// return vector
	d_vector val_out(K);
	//
	// compute the sparse Hessian
	ran_like_fun_.SparseHessian(
		both,
		w,
		not_used,
		ran_hes_.row,
		ran_hes_.col,
		val_out,
		ran_hes_.work
	);
	//
	double eps = 100. * std::numeric_limits<double>::epsilon();
	bool ok = true;
	for(size_t k = 0; k < K; k++)
	{	ok &= CppAD::NearEqual(
			Value(Var2Par(a1_val_out[k])), val_out[k], eps, eps
		);
	}
	if( ! ok )
	{	const std::string error_message =
		"init_ran_hes: a ran_likelihood_hes Hessian value "
		"does not agree with AD value.";
		fatal_error(error_message);
	}
	return;
}
