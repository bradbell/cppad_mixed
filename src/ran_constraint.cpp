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
------------------------------------------------------------------------------
$begin init_ran_con_mat$$
$spell
$$

$section Initialize Random Constraints Matrix$$

$head Syntax$$
$icode%mixed_object%.init_ran_con_mat(%n_random%, %A_info%)
%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is the number of
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$.
It is assumed that $icode%n_random% > 0%$$.

$head A_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %A_info%
%$$
It is a
$cref/sparse matrix/sparse_mat_info/Sparse Matrix/$$ representation
of the
$cref/random constraint matrix
	/cppad_mixed
	/Notation
	/Random Constraint Matrix, A
/$$
$latex A$$.

$subhead Restriction$$
We we assume that $icode%A_info.val[%k%] != 0%$$ for
$icode%k% = 0 , ... , %K%-1%$$ where
$icode%K% = %A_info%.row.size()%$$.

$head ran_con_mat_$$
The input value of the member variable
$codei%
	CppAD::mixed::cholesky::eigen_sparse ran_con_mat_
%$$
does not matter.
Upon return it is an $code Eigen$$ sparse matrix representation
of the random constraint matrix $latex A$$.

$end
*/

void cppad_mixed::init_ran_con_mat(
	size_t                               n_random ,
	const CppAD::mixed::sparse_mat_info& A_info   )
{
	// number of possibly non-zero entries
	size_t K = A_info.row.size();
	//
	assert( n_random > 0 );
	assert( A_info.col.size() == K );
	assert( A_info.val.size() == K );

	// determine maximum row and column index
	size_t max_row_p1 = 0;
	size_t max_col_p1 = 0;
	for(size_t k = 0; k < K; k++)
	{	assert( A_info.val[k] != 0.0 );
		max_row_p1 = std::max(max_row_p1, A_info.row[k] + 1);
		max_col_p1 = std::max(max_col,    A_info.col[k] + 1);
	}
	assert( max_col_p1 <= n_random );

	// size of the matrix A
	nrow = max_row_p1;
	ncol = n_random;
	ran_con_mat_.resize(nrow, ncol);
	// we know number of non-zeros
	ran_con_mat_.reserve(K);

	// put the values in A
	for(size_t k = 0; k < K; k++)
		ran_con_mat_.insert(A_info.row[k], A_info.col[k]) = A_info.val[k];

	return;
}

/*
-------------------------------------------------------------------------------
$begin ran_con_eval$$
$spell
$$

$section Evaluate the Random Constraint Function$$

$head Syntax$$
$icode%vec% = mixed_object%.ran_con_eval(%random_vec%)
%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head random_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_vec%
%$$
It specifies the value of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$ at which $icode%A%*u%$$ is evaluated.

$head vec$$
The return value has prototype
$codei%
	CppAD::vector<double> %vec%
%$$
The $icode random_vec$$ is the
$cref/optimal random effects
	/theory
	/Optimal Random Effects, u^(theta)
/$$
$latex \hat{u} ( \theta )$$,
$icode vec$$ is the value of the
$cref/random constraint Function
	/cppad_mixed
	/Notation
	/Random Constraint Function, A*u
/$$
$latex A \hat{u) ( \theta ) $$.

$children%
	example/private/ran_con_eval_xam.cpp
%$$
$head Example$$
The file $cref ran_con_eval_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
void CppAD::vector<double> cppad_mixed::ran_con_eval(
	const d_vector& random_vec )
{	assert( random_vec.size() == n_random_ );
	//
	// copy random_vec to an eigen column matrix
	typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_dense;
	eigen_dense u( n_random_ , 1);
	for(size_t j = 0; j < n_random_; j++)
		u[j] = random_vec[j];
	//
	// multiply by the constraint matrix
	eigen_dense Au = ran_con_mat_ * u;
	//
	// copy the return value
	d_vector vec( ran_con_mat_.rows() );
	for(size_t i = 0; i < vec.size(); i++)
		vec[i] = Au[i];
	//
	return vec;
}
/*
-------------------------------------------------------------------------------
$begin ran_con_jac$$
$spell
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
It is a
$cref/sparse matrix/sparse_mat_info/Sparse Matrix/$$ representation
of the Jacobian of the
$cref/random constraint function
	/cppad_mixed
	/Notation
	/Random Constraint Function, A*u
/$$; i.e.,
$latex A \hat{u} ( \theta )$$; i.e., the matrix
$latex \[
	\partial_\theta [ A \hat{u} ( \theta ) ]
	=
	A \hat{u}_\theta ( \theta )
	=
	- A
	f_{u,u} \left[ \theta , \hat{u} ( \theta ) \right]^{-1}
	f_{u,\theta} \left[ \theta , \hat{u} ( \theta )  \right]
\] $
see
$cref/derivative of optimal random effects
	/theory
	/Derivative of Optimal Random Effects
/$$.

$children%
	example/private/ran_con_eval_xam.cpp
%$$
$head Example$$
The file $cref ran_con_eval_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
void CppAD::vector<double> cppad_mixed::ran_con_eval(
	const d_vector&          fixed_vec  ,
	const d_vector&          random_vec ,
	CppAD::vector<size_t>&   row        ,
	CppAD::vector<size_t>&   col        ,
	d_vector&                val        )
{

}

