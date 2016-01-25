// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin ran_con_jac$$
$spell
	Jacobian
	jac
	vec
	cppad
	const
	CppAD
$$

$section Jacobian of the Random Constraint Function$$

$head Syntax$$
$icode%mixed_object%.ran_con_jac(
	%fixed_vec%, %random_vec%, %jac_info%
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
vector $latex \theta$$.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$,
This must be the
$cref/optimal random effects
	/theory
	/Optimal Random Effects, u^(theta)
/$$
$latex \hat{u} ( \theta )$$.

$head jac_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %jac_info%
%$$
The input size of the vectors in $icode jac_info$$
is either zero, or the vectors are the same as the return value
for a previous call to $code ran_con_jac$$.

$subhead Sparsity Patter$$
In the case where the input size of the vectors in $icode jac_info$$
is zero, only the
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$ is computed.
Upon return, the size of $icode%jac_info%.val%$$ is equal to the
size of $icode%jac_info%.row%$$, but the elements of
$icode%jac_info%.val%$$ are not specified.

$subhead Sparse Matrix$$
If the input size of the vectors in $icode jac_info$$
are non-zero,
upon return $icode jac_info$$ is a
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
representation of the Jacobian of the
$cref/random constraint function
	/cppad_mixed
	/Notation
	/Random Constraint Function, A*u^(theta)
/$$; i.e.,
$latex \[
	\partial_\theta [ A \hat{u} ( \theta ) ]
	=
	A \hat{u}_\theta ( \theta )
	=
	- A
	f_{u,u} \left[ \theta , \hat{u} ( \theta ) \right]^{-1}
	f_{u,\theta} \left[ \theta , \hat{u} ( \theta )  \right]
\] $$
see
$cref/derivative of optimal random effects
	/theory
	/Derivative of Optimal Random Effects
/$$.

$subhead Column Major Order$$
The results in $icode jac_info$$ are in
$cref/column major order/sparse_mat_info/Notation/Column Major Order/$$.

$children%
	example/private/ran_con_jac_xam.cpp
%$$

$head Example$$
The file $cref ran_con_jac_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::ran_con_jac(
	const d_vector&                fixed_vec  ,
	const d_vector&                random_vec ,
	CppAD::mixed::sparse_mat_info& jac_info   )
{	assert( fixed_vec.size()  == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	//
	size_t L = jac_info.row.size();
	assert( jac_info.col.size() == L );
	assert( jac_info.val.size() == L );

	// declare eigen matrix types
	typedef typename CppAD::mixed::choleig::eigen_sparse eigen_sparse;
	typedef typename eigen_sparse::InnerIterator          column_itr;

	// packed version of fixed and random effects
	d_vector both(n_fixed_ + n_random_);
	pack(fixed_vec, random_vec, both);

	// Compute the Hessian cross terms f_{u,theta} ( theta , u )
	CppAD::vector< std::set<size_t> > not_used;
	size_t K = hes_cross_.row.size();
	d_vector val_out(K), w(1);
	w[0] = 1.0;
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
	// initialize matrix product A * uhat_theta( \theta )
	assert( n_ran_con_ == size_t ( ran_con_mat_.rows() ) );
	eigen_sparse Au_theta( n_ran_con_, n_fixed_ );

	// Loop over fixed effects
	size_t ell = 0;
	CppAD::vector<size_t> row_b, row_x;
	CppAD::vector<double> val_b, val_x;
	for(size_t j = 0; j < n_fixed_; j++)
	{	// b = j-th column of - f_{u, theta} (theta, u)
		row_b.resize(0);
		val_b.resize(0);
		while( col <= j )
		{	assert( col == j );
			// Must put a non-zero value on left hand size
			// because Eigens Cholesky solve seems to drop zeros.
			if( L == 0 )
			{	row_b.push_back(row);
				val_b.push_back(1.0);
			}
			else
			{	row_b.push_back(row);
				val_b.push_back( - val_out[k] );
			}
			//
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
		// x = j-th column of - f_{u,u}(theta, u)^{-1} f_{u,theta}(theta, u)
		row_x.resize(0);
		val_x.resize(0);
		chol_ran_hes_.solve(row_b, val_b, row_x, val_x);
		//
		// convert to an eigen column vector
		eigen_sparse x(n_random_, 1);
		for(size_t m = 0; m < row_x.size(); m++)
			x.insert( row_x[m], 0) = val_x[m];
		//
		// multipliy A times j-th column of uhat_theta
		eigen_sparse jac_j = ran_con_mat_ * x;
		// convert to sparse_mat_info format
		for(column_itr itr(jac_j, 0); itr; ++itr)
		{	size_t i = itr.row();
			if( L == 0 )
			{	jac_info.row.push_back(i);
				jac_info.col.push_back(j);
				jac_info.val.push_back(0.0);
			}
			else
			{	// check for case where eigen has dropped a zero elements
				while( jac_info.row[ell] < i )
					jac_info.val[ell++] = 0.0;
				//
				assert( jac_info.row[ell] == i );
				assert( jac_info.col[ell] == j );
				jac_info.val[ell] = itr.value();
				ell++;
			}
		}
	}
	return;
}
