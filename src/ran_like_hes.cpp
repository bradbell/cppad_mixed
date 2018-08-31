// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
$begin ran_like_hes$$
$spell
	hes_rcv
	jac
	cppad
	vec
	const
	xam
	init
$$

$section Hessian of Random Likelihood w.r.t. Random Effects$$

$head Syntax$$
$icode%hes_rcv% = %mixed_object%.ran_like_hes(%fixed_vec%, %random_vec%)%$$

$head Prototype$$
$srcfile%src/ran_like_hes.cpp
    %0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Assumptions$$
The member variable
$cref/init_ran_jac_done_/init_ran_jac/init_ran_jac_done_/$$ is true.

$head Purpose$$
This routine computes the Hessian of the random likelihood
$cref/f(theta, u)/theory/Random Likelihood, f(theta, u)/$$
with respect to the random effects vector $latex u$$; i.e.
$latex \[
	f_{uu} ( \theta, u )
\] $$

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_vec$$
This argument specifies the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.

$head random_vec$$
This argument specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.

$head hes_rcv$$
The return value
is the Hessian $latex f_{uu} ( \theta , u )$$.
The indices in this matrix are just with respect to the random effects;
i.e., the row and column indices are between zero and  $code n_random_$$.

$children%
	example/private/ran_like_hes.cpp
%$$
$head Example$$
The file $cref ran_like_hes.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.
----------------------------------------------------------------------------
$end
*/
// BEGIN_PROTOTYPE
CppAD::mixed::a1_sparse_rcv cppad_mixed::ran_like_hes(
	const a1_vector& fixed_vec   ,
	const a1_vector& random_vec  )
// END_PROTOTYPE
{	assert( init_ran_jac_done_ );
	assert( init_ran_hes_done_ );
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_vec.size() );
	//
	// ran_hes_uu_rc
	const sparse_rc& ran_hes_uu_rc( a1_ldlt_ran_hes_.pattern() );
	//
	// row, col
	const s_vector&  row( ran_hes_uu_rc.row() );
	const s_vector&  col( ran_hes_uu_rc.col() );
	//
	// nnz
	size_t nnz = ran_hes_uu_rc.nnz();
	//
	// The user has not defined ran_likelihood_hes, so use AD to calcuate
	// the Hessian of the random likelihood w.r.t the random effects.
	//
	// ran_hes_mix_rc
	sparse_rc ran_hes_mix_rc(n_random_, n_fixed_ + n_random_, nnz);
	for(size_t k = 0; k < nnz; k++)
	{	assert( row[k] <= n_random_ );
		assert( col[k] <= n_random_ );
		assert( col[k] <= row[k] );
		//
		// row relative to random effects, column relative to both
		ran_hes_mix_rc.set(k, row[k], col[k] + n_fixed_);
	}
	// theta_u
	a1_vector theta_u(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, theta_u);
	//
	// ran_hes_mix_rcv
	a1_sparse_rcv ran_hes_mix_rcv( ran_hes_mix_rc );
	ran_jac_a1fun_.subgraph_jac_rev(theta_u, ran_hes_mix_rcv);
	ran_jac_a1fun_.clear_subgraph();
	//
	// ran_hes_uu_rcv
	a1_sparse_rcv ran_hes_uu_rcv( ran_hes_uu_rc );
	for(size_t k = 0; k < nnz; ++k)
		ran_hes_uu_rcv.set(k, ran_hes_mix_rcv.val()[k] );
	//
	return ran_hes_uu_rcv;
}
