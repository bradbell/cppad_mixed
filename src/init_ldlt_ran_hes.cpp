// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-17 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin init_ldlt_ran_hes$$
$spell
	rcv
	ldlt_ran_hes
	init
	cppad
	const
	CppAD
	cholesky
$$

$section Initialize Cholesky Factor of Hessian of Random Likelihood$$

$head Syntax$$
$icode%mixed_object%.init_ldlt_ran_hes()%$$

$head Private$$
This $code cppad_mixed$$ member function is $cref private$$.

$head Assumptions$$
The member variable
$cref/init_ran_hes_done_/init_ran_hes/init_ran_hes_done_/$$ is true.

$head init_ldlt_ran_hes_done_$$
The input value of this member variable must be false.
Upon return it is true.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head ldlt_ran_hes_$$
The member variable
$codei%
	CPPAD_MIXED_LDLT_CLASS ldlt_ran_hes_
%$$
must not been previously initialized
$codei%
	size_t ldlt_ran_hes_
%$$
Upon return, the function
$codei%
	ldlt_ran_hes_.init(%hes_info%)
%$$
has been called with $icode hes_info$$
equal to the sparsity pattern for the
Hessian of the random likelihood with respect to the random effects;
see $cref ldlt_eigen_init$$.
This sparsity information is relative to just the random effects;
i.e., the row and column indices are relative to the vector $latex u )$$.
For each row and column indices in $icode hes_info$$,
the corresponding row and column index in $code ran_hes_rcv_$$ is
$code n_fixed_$$ greater.


$end
*/
# include <cppad/mixed/cppad_mixed.hpp>

void cppad_mixed::init_ldlt_ran_hes(void)
{	assert( ! init_ldlt_ran_hes_done_ );
	assert( init_ran_hes_done_ );
	//
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
	ldlt_ran_hes_.init(hes_info);
	//
	init_ldlt_ran_hes_done_ = true;
	return;
}
