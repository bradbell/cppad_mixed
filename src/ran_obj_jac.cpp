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
$begin ran_obj_jac$$
$spell
	chol
	jac
	CppAD
	ran_obj
	cppad
	hes
	vec
	const
	Cpp
	xam
	cholesky
$$

$section Derivative of Random Objective$$

$head Syntax$$
$icode%mixed_object%.ran_obj_jac(%fixed_vec%, %random_vec%, %r_fixed%)%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Purpose$$
This routine computes the
$cref/derivative of the random objective
	/theory/
	Derivative of Random Objective
/$$.

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
	example/private/ran_obj_jac_xam.cpp
%$$
$head Example$$
The file $cref ran_obj_jac_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------
void cppad_mixed::ran_obj_jac(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec ,
	d_vector&       r_fixed    )
{
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( r_fixed.size() == n_fixed_ );

	// number of non-zeros in Hessian
	size_t K = ran_hes_.row.size();
	assert( K == ran_hes_.col.size() );

	//
	// Compute derivative of logdet of f_{u,u} ( theta , u )
	d_vector logdet_fix(n_fixed_), logdet_ran(n_random_);
	logdet_jac(fixed_vec, random_vec, logdet_fix, logdet_ran);

	// packed version of fixed and random effects
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);

	//
	// Compute derivative of f(theta , u ) w.r.t theta and u
	d_vector w(1), f_both(n_fixed_ + n_random_);
	w[0] = 1.0;
	ran_like_fun_.Forward(0, both);
	f_both = ran_like_fun_.Reverse(1, w);
	d_vector f_fixed(n_fixed_), f_random(n_random_);
	unpack(f_fixed, f_random, f_both);
	//
	// Compute the Hessian cross terms f_{u,theta} ( theta , u )
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

	// Use the column major order specification for
	// (hes_cross_.row, hes_cross_.col)
	size_t k = 0;
	size_t row = n_random_;
	size_t col = n_fixed_;
	if( k < K )
	{	assert( hes_cross_.row[k] >= n_fixed_ );
		row = hes_cross_.row[k] - n_fixed_;
		col = hes_cross_.col[k];
		assert( row < n_random_ );
		assert( col < n_fixed_ );
	}

	// Loop over fixed effects and compute r_fixed one component at a time
	CppAD::vector<bool> found(n_random_);
	CppAD::vector<size_t> row_solve;
	CppAD::vector<double> val_b, val_x;
	for(size_t j = 0; j < n_fixed_; j++)
	{	// vectors for this column
		row_solve.resize(0);
		val_b.resize(0);
		val_x.resize(0);
		//
		// increment through all row indices
		size_t i = 0;

		// b = j-th column of - f_{u, theta} (theta, u)
		while( col <= j )
		{	// make sure we are solving for all rows where
			// derivative of logdet w.r.t random effects is non-zero.
			while( i < row )
			{	if( logdet_ran[i] != 0.0 )
				{	row_solve.push_back(i);
					val_b.push_back(0.0);
				}
				i++;
			}
			if( col == j )
			{	row_solve.push_back(row);
				val_b.push_back( - val_out[k] );
				found[row] = true;
				assert( i == row );
				i++;
			}
			k++;
			if( k < K )
			{	assert( hes_cross_.row[k] >= n_fixed_ );
				row = hes_cross_.row[k] - n_fixed_;
				col = hes_cross_.col[k];
				assert( row < n_random_ );
				assert( col < n_fixed_ );
			}
			else
			{	row = n_random_;
				col = n_fixed_;
			}
		}
		assert( col > j );
		// make sure we are solving for all rows where
		// derivative of logdet w.r.t random effects is non-zero.
		while( i < n_random_ )
		{	if( logdet_ran[i] != 0.0 )
			{	row_solve.push_back(i);
				val_b.push_back(0.0);
			}
			i++;
		}

		// x = j-th column of - f_{u,u}(theta, u)^{-1} f_{u,theta}(theta, u)
		val_x.resize( row_solve.size() );
		chol_ran_hes_.solve_H(row_solve, val_b, val_x);
		//
		// parial w.r.t fixed effects contribution to total derivative
		r_fixed[j] =  f_fixed[j] + 0.5 * logdet_fix[j];
		//
		// compute effect of uhat_{theta(j)} (theta) on the total derivative
		for(size_t ell = 0; ell < row_solve.size(); ell++)
		{	// random effect index
			size_t i = row_solve[ell];
			// partial of optimal random effect for this (i, j)
			double ui_thetaj = val_x[ell];
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
