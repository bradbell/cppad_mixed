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
# include <cppad/mixed/exception.hpp>
/*
$begin update_factor$$
$spell
	Taylor
	vec
	Cholesky
	cppad
	const
	CppAD
	ldlt
	hes
$$

$section Update the Factorization of Hessian w.r.t. Random Effects$$

$head Syntax$$
$icode%mixed_object%.update_factor(%fixed_vec%, %random_vec%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Purpose$$
This routine updates the Cholesky factorization of the Hessian
$latex f_{u,u} ( \theta , u )^{-1} $$
so that it corresponds to the current value of
$latex \theta$$ and $latex u$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
and is the value of fixed effects $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
and is the value of fixed effects $latex u$$.
If this factorization is used for computing derivatives of the
$cref/random objective/theory/Objective/Random Objective, r(theta)/$$,
this should be the
$cref/optimal random effects/theory/Optimal Random Effects, u^(theta)/$$.

$head ran_hes_fun_$$
The $cref/ran_hes_fun_/private/ran_hes_fun_/$$ member variable
will hold the first order Taylor coefficient corresponding
to the specified fixed and random effects; i.e.,
$codei%
	ran_hes_fun_.Forward(0, %both%)
%$$
has been called where $icode both$$ is a packed version
of the fixed and random effects.

$head ldlt_ran_hes_$$
On input, the member variable
$codei%
	CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_
%$$
has been
$cref/initialized/ldlt_eigen_init/$$
using the sparsity pattern for the Hessian.
Upon return, $code ldlt_ran_hes_$$ contains the updated
$cref/factorization/ldlt_eigen_update/$$
corresponding to the specified values for the fixed
and random effects.

$children%
	example/private/update_factor.cpp
%$$
$head Example$$
The file $cref update_factor.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------
void cppad_mixed::update_factor(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec )
{	assert( init_ran_hes_done_ );
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	//
	// set the sparsity pattern corresponding the Hessian
	CppAD::mixed::sparse_mat_info hes_info;
	size_t K = ran_hes_rcv_.nnz();
	hes_info.row.resize(K);
	hes_info.col.resize(K);
	for(size_t k = 0; k < K; k++)
	{	size_t r = ran_hes_rcv_.row()[k];
		size_t c = ran_hes_rcv_.col()[k];
		assert( n_fixed_ <= r && r < n_fixed_ + n_random_ );
		assert( n_fixed_ <= c && c < n_fixed_ + n_random_ );
		hes_info.row[k] = r - n_fixed_;
		hes_info.col[k] = c - n_fixed_;
	}
	//
	// pack fixed and random effects into one vector
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);
	//
	// set the value vector in the sparse matrix information
	hes_info.val = ran_hes_fun_.Forward(0, both);
	if( CppAD::hasnan( hes_info.val ) ) throw CppAD::mixed::exception(
		"update_factor", "result has nan"
	);
	//
	// update the LDLT factor
	bool ok = ldlt_ran_hes_.update(hes_info);
	if( ! ok )
	{
# if CPPAD_MIXED_LOG_FATAL_ERREOR
		CppAD::mixed::exception e(
			"update_factor", "Hessian w.r.t. random effects is singular"
		);
		throw(e);
# else
		// convert fatal error to an assert (for use in a debugger)
		assert(false);
# endif
	}
	return;
}

