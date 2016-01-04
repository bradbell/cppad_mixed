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
# include <cppad/mixed/chol_hes_ran.hpp>
/*
$begin ran_obj_grad$$
$spell
	CppAD
	ran_obj
	cppad
	hes
	vec
	const
	Cpp
	xam
$$

$section Derivative of Random Objective$$

$head Syntax$$
$icode%mixed_object%.ran_obj_grad(%fixed_vec%, %random_vec%, %r_fixed%)%$$

$head Purpose$$
This routine computes the
$cref/derivative of the random objective
	/theory/
	Derivative of Random Objective
/$$.

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
and is the value of fixed effects $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
and is the value of fixed effects $latex u$$.
It is assumed that these random effects are optimal for the specified
fixed effects and hence $latex f_u ( \theta , u ) = 0$$.

$head r_fixed$$
This argument has prototype
$codei%
	CppAD::vector<double>& %r_fixed%
%$$
If the input size must be equal to $code n_fixed_$$.
Upon return, it contains the value of the derivative w.r.t
the fixed effects; i.e. $latex r_\theta ( \theta )$$.

$children%
	example/private/ran_obj_grad_xam.cpp
%$$
$head Example$$
The file $cref ran_obj_grad_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------
void cppad_mixed::ran_obj_grad(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec ,
	d_vector&       r_fixed    )
{
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( r_fixed.size() == n_fixed_ );

	// declare eigen matrix types
	typedef Eigen::SparseMatrix<double>                sparse_matrix;
	typedef Eigen::SparseMatrix<double>::InnerIterator inner_itr;

	// number of non-zeros in Hessian
	size_t K = hes_ran_.row.size();
	assert( K == hes_ran_.col.size() );

	// compute an LDL^T Cholesky factorization of f_{u,u}(theta, u)
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);
	CppAD::mixed::factorize_chol_hes_ran(
		n_fixed_, n_random_, hes_ran_.row, hes_ran_.col, both, hes_ran_fun_
	);

	//
	// Compute derivative of logdet of f_{u,u} ( theta , u )
	d_vector logdet_fix(n_fixed_), logdet_ran(n_random_);
	logdet_jac(fixed_vec, random_vec, logdet_fix, logdet_ran);
	//
	// Compute derivative of f(theta , u ) w.r.t theta and u
	d_vector w(1), f_both(n_fixed_ + n_random_);
	w[0] = 1.0;
	ran_like_fun_.Forward(0, both);
	f_both = ran_like_fun_.Reverse(1, w);
	d_vector f_fixed(n_fixed_), f_random(n_random_);
	unpack(f_fixed, f_random, f_both);
	//
	// Compute the Hessian cross terms f_{u theta}^{(2)} ( theta , u )
	// 2DO: another ran_like_fun_.Forward(0, both) is done by SparseHessian
	CppAD::vector< std::set<size_t> > not_used;
	K = hes_cross_.row.size();
	CppAD::vector<double> val_out(K);
	ran_like_fun_.SparseHessian(
		both,
		w,
		not_used,
		hes_cross_.row,
		hes_cross_.col,
		val_out,
		hes_cross_.work
	);
	// initialize index in hes_cross_.row, hes_cross_.col
	size_t k = 0;
	size_t col = n_fixed_;
	if( k < K )
		col = hes_cross_.col[k];
	assert( col <= n_fixed_ );
	//
	// Loop over fixed effects and compute r_fixed one component at a time
	for(size_t j = 0; j < n_fixed_; j++)
	{	// parial w.r.t fixed effects contribution to total derivative
		r_fixed[j] =  f_fixed[j] + 0.5 * logdet_fix[j];
		//
		// set b_i = - f_{u theta}^{(2)} (theta , u) ]_{i,j}
		sparse_matrix b(n_random_, 1);
		while( col < j )
		{	k++;
			if( k >= K )
				col = n_fixed_;
			else
			{	col = hes_cross_.col[k];
				assert( col < n_fixed_ );
			}
		}
		while( col == j )
		{	assert( hes_cross_.row[k] >= n_fixed_ );
			size_t row = hes_cross_.row[k] - n_fixed_;
			assert( row < n_random_ );
			b.insert(row, 0) = - val_out[k];
			k++;
			if( k >= K )
				col = n_fixed_;
			else
			{	col = hes_cross_.col[k];
				assert( col < n_fixed_ );
			}
		}
		//
		// compute the partial of uhat(theta) w.r.t theta[j]
		sparse_matrix x = CppAD::mixed::chol_hes_ran_.solve(b);
		//
		// compute effect on the total derivative
		assert( x.outerSize() == 1 );
		for(inner_itr itr(x, 0); itr; ++itr)
		{	// random effect index
			size_t i = itr.row();
			// partial of optimal random effect for this (i, j)
			double ui_thetaj = itr.value();
			// partial of random part of objective w.r.t this random effect
			// Note f_random = 0 because u is optimal for the fixed effects
			// double h_ui = f_random[i] + 0.5 * logdet_ran[i];
			double h_ui = 0.5 * logdet_ran[i];
			// contribution to total derivative
			r_fixed[j] += h_ui * ui_thetaj;
		}
	}
	//
	return;
}

