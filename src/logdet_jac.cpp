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
$begin logdet_jac$$
$spell
	chol
	Jacobian
	jac
	CppAD
	cppad
	hes
	vec
	const
	Cpp
	xam
	logdet
	cholesky
$$

$section Jacobian of Log Determinant of Hessian w.r.t. Random Effects$$

$head Syntax$$
$icode%mixed_object%.logdet_jac(
	%fixed_vec%, %random_vec%, %logdet_fix%, %logdet_ran%
)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Purpose$$
This routine computes the Jacobian of the log determinant
of the Hessian of the random likelihood
$cref/f(theta, u)
	/theory/
	Random Likelihood, f(theta, u)
/$$
with respect to the random effects vector $latex u$$.
To be specific, it computes both
$latex \[
	\partial_\theta \log \det [ f_{u,u} ( \theta, u ) ]
	\; \R{and} \;
	\partial_u \log \det [ f_{u,u} ( \theta, u ) ]
\] $$

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head chol_ran_hes_$$
It is assumed that the member variable
$codei%
	CPPAD_MIXED_LDLT chol_ran_hes_
%$$
was updated using $cref update_factor$$ for the specified values of the
fixed and random effects.

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
$codei%ran_hes_.col[%k%] - n_fixed_%$$ of the Hessian.

$head logdet_ran$$
This argument has prototype
$codei%
	CppAD::vector<double>& %logdet_ran%
%$$
Its input size must be equal to $code n_random_$$.
Upon return, it contains the value of the derivative w.r.t
the random effects.

$children%
	example/private/logdet_jac_xam.cpp
%$$
$head Example$$
The file $cref logdet_jac_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------
void cppad_mixed::logdet_jac(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec ,
	d_vector&       logdet_fix ,
	d_vector&       logdet_ran )
{	assert( init_ran_hes_done_ );

	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( logdet_fix.size() == n_fixed_ );
	assert( logdet_ran.size() == n_random_ );

	// number of non-zeros in Hessian
	size_t K = ran_hes_.row.size();
	assert( K == ran_hes_.col.size() );

	// Compute derivative of sum_k w_k hessian_k
	// where w_k is f_{u,u} (theta, u)^{-1} at (row[k], col[k]).
	d_vector w(K);
	for(size_t k = 0; k < K; k++)
		w[k] = 0.0;

	// starting index in (ran_hes_.row, ran_hes_.col)
	size_t k  = 0;
	size_t col = n_random_;
	size_t row = n_random_;
	if( k < K )
	{	assert( ran_hes_.row[k] >= n_fixed_ );
		assert( ran_hes_.col[k] >= n_fixed_ );
		row = ran_hes_.row[k] - n_fixed_;
		col = ran_hes_.col[k] - n_fixed_;
		assert( row < n_random_ );
		assert( col < n_random_ );
	}
	CppAD::vector<size_t> row_solve;
	d_vector val_in, val_out;
	for(size_t j = 0; j < n_random_; j++)
	{	// vectors for this column
		row_solve.resize(0);
		val_in.resize(0);
		val_out.resize(0);

		// only need for rows where f_{u,u} (theta_u) is possibly not zero
		size_t k_start    = K;
		bool   row_j_zero = true;
		while( col <= j )
		{	if( col == j )
			{	// this row of f_{u,u} (theta, u) is possibly non-zero
				row_solve.push_back( row );
				// val_in needs to be j-th column of identity matrix
				// so solution is j-th column of inverse
				if( row == j )
				{	val_in.push_back(1.0);
					row_j_zero = false;
				}
				else
					val_in.push_back(0.0);
				//
				// index in hes_ran_ where j-th column starts
				if( k_start == K )
					k_start = k;
			}
			k++;
			if( k < K )
			{	assert( ran_hes_.row[k] >= n_fixed_ );
				assert( ran_hes_.col[k] >= n_fixed_ );
				row = ran_hes_.row[k] - n_fixed_;
				col = ran_hes_.col[k] - n_fixed_;
				assert( row < n_random_ );
				assert( col < n_random_ );
			}
			else
				row = col = n_random_;
		}
		assert( col > j );
		//
		// Cannot compute cholesky factor if f_{u,u} (theta, u) is zero
		assert( ! row_j_zero );
		if( row_j_zero )
		{	row_solve.push_back(j);
			val_in.push_back(1.0);
		}
		// if k_start == K, we do not need any components of the inverse
		if( k_start < K )
		{	val_out.resize( row_solve.size() );
			//
			chol_ran_hes_.solve_H(row_solve, val_in, val_out);
			//
			size_t nrow = row_solve.size();
			if( row_j_zero )
				nrow--;
			for(size_t ell = 0; ell < nrow; ell++)
			{	size_t k    = k_start + ell;
				// note off diagonal elements need to be counted twice
				// becasue only computing for lower triangle
				w[k] = val_out[ell];
				if( row_solve[ell] != j )
					w[k] = 2.0 * w[k];
			}
		}
	}
	d_vector dw(n_fixed_ + n_random_);
	dw = ran_hes_fun_.Reverse(1, w);
	//
	// split out fixed and random parts of the derivative
	unpack(logdet_fix, logdet_ran, dw);
	//
	return;
}

