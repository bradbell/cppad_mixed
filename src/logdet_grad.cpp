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
# include <dismod_at/chol_hes_ran.hpp>
/*
$begin logdet_grad$$
$spell
	cppad
	hes
	vec
	const
	Cpp
	xam
	logdet
$$

$section Derivative of Log Determinant of Hessian w.r.t. Random Effects$$

$head Syntax$$
$icode%mixed_object%.logdet_grad(
	%fixed_vec%, %random_vec%, %logdet_fix%, %logdet_ran%)%$$

$head Purpose$$
This routine computes the total derivative of the log determinant
of the Hessian of the random likelihood
$cref/f(theta, u)
	/cppad_mixed_theory/
	Random Likelihood, f(theta, u)
/$$
with respect to the random effects vector $latex u$$; i.e.
it computes both
$latex \[
	\partial_\theta \log \det [ f_{uu}^{(2)} ( \theta, u ) ]
	\; \R{and} \;
	\partial_u \log \det [ f_{uu}^{(2)} ( \theta, u ) ]
\] $$

$head mixed_object$$
We use $cref/mixed_object/cppad_mixed_derived_ctor/mixed_object/$$
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

$head logdet_fix$$
This argument has prototype
$codei%
	CppAD::vector<double>& %logdet_fix%
%$$
Its input size must be equal to $code n_fixed_$$.
Upon return, it contains the value of the derivative w.r.t
the fixed effects.
$codei%hes_ran_.col[%k%] - n_fixed_%$$ of the Hessian.

$head logdet_ran$$
This argument has prototype
$codei%
	CppAD::vector<double>& %logdet_ran%
%$$
Its input size must be equal to $code n_random_$$.
Upon return, it contains the value of the derivative w.r.t
the random effects.

$children%
	example/devel/cppad_mixed/private/logdet_grad_xam.cpp
%$$
$head Example$$
The file $cref logdet_grad_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
namespace dismod_at { // BEGIN_DISMOD_AT_NAMESPACE
// ----------------------------------------------------------------------------
void cppad_mixed::logdet_grad(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec ,
	d_vector&       logdet_fix ,
	d_vector&       logdet_ran )
{	assert( init_hes_ran_done_ );

	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( logdet_fix.size() == n_fixed_ );
	assert( logdet_ran.size() == n_random_ );

	// declare eigen matrix types
	typedef Eigen::SparseMatrix<double>                sparse_matrix;
	typedef Eigen::SparseMatrix<double>::InnerIterator inner_itr;

	// number of non-zeros in Hessian
	size_t K = hes_ran_.row.size();
	assert( K == hes_ran_.col.size() );

	// compute an LDL^T Cholesky factorization of f_{uu}^{(2)}(theta, u)
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);
	factorize_chol_hes_ran(
		n_fixed_, n_random_, hes_ran_.row, hes_ran_.col, both, hes_ran_fun_
	);

	// Compute derivative of sum_k w_k hessian_k
	// where w_k is the inverse of the Hessian at (row[k], col[k]).
	// use order speicifcation for (hes_ran_.row, hes_ran_.col)
	d_vector w(K);
	for(size_t k = 0; k < K; k++)
		w[k] = 0.0;
	size_t k  = 0;
	size_t col = n_random_;
	size_t row = n_random_;
	if( k < K )
	{	row = hes_ran_.row[k] - n_fixed_;
		col = hes_ran_.col[k] - n_fixed_;
	}
	for(size_t j = 0; j < n_random_; j++)
	{	while( col < j )
		{	k++;
			if( k < K )
			{	row = hes_ran_.row[k] - n_fixed_;
				col = hes_ran_.col[k] - n_fixed_;
			}
			else
				row = col = n_random_;
		}
		sparse_matrix b(n_random_, 1);
		b.insert(j, 0) = 1.0;
		sparse_matrix x = chol_hes_ran_.solve(b);
		assert( x.outerSize() == 1 );
		for(inner_itr itr(x, 0); itr; ++itr)
		{	size_t i    = itr.row();
			if( col == j )
			{	while( row < i )
				{	k++;
					if( k < K )
					{	row = hes_ran_.row[k] - n_fixed_;
						col = hes_ran_.col[k] - n_fixed_;
					}
					else
						row = col = n_random_;
				}
			}
			if( (row == i) && col == j)
			{	// note off diagonal elements need to be counted twice
				// becasue only computing for lower triangle
				w[k] = itr.value();
				if( hes_ran_.row[k] != hes_ran_.col[k] )
					w[k] = 2.0 * w[k];
			}
		}
	}
	d_vector dw(n_fixed_ + n_random_);
	dw = hes_ran_fun_.Reverse(1, w);
	//
	// split out fixed and random parts of the derivative
	unpack(logdet_fix, logdet_ran, dw);
	//
	return;
}

} // END_DISMOD_AT_NAMESPACE
