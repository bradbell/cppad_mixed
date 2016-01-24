// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */

# include <Eigen/Sparse>
# include <cppad/mixed/cholesky.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// constructor
cholesky::cholesky(size_t n_random)
: n_random_(n_random)
{	ptr_ = new eigen_cholesky; }

// destructor
cholesky::~cholesky(void)
{	delete ptr_; }

/*
------------------------------------------------------------------------------
$begin cholesky_init$$
$spell
	Cholesky
	chol_ran_hes
	CppAD
	const
	init
$$

$section Initialize Cholesky Factor for a Specific Sparsity Pattern$$

$head Syntax$$
$icode%chol_ran_hes%.init(%hes_info%)
%$$

$head Private$$
The $code cholesky$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head chol_ran_hes$$
This object has prototype
$codei%
	CppAD::mixed::cholesky %chol_ran_hes%
%$$

$head hes_info$$
This argument had prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It is a
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$ for the
matrix that we will compute the Cholesky factor of.
It is in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$end
*/
void cholesky::init( const CppAD::mixed::sparse_mat_info& hes_info )
{	double not_used = 1.0;

	Eigen::SparseMatrix<double> hessian_pattern(n_random_, n_random_);
	assert( hes_info.row.size() == hes_info.col.size() );
	for(size_t k = 0; k < hes_info.row.size(); k++)
	{	size_t r = hes_info.row[k];
		size_t c = hes_info.col[k];
		assert( r < n_random_ );
		assert( c < n_random_ );
		hessian_pattern.insert(r, c) = not_used;
	}
	// analyze the pattern for an LDL^T Cholesky factorization of
	// f_{u,u}(theta, u)
	ptr_->analyzePattern(hessian_pattern);
}
/*
------------------------------------------------------------------------------
$begin cholesky_update$$
$spell
	Cholesky
	xam
	const
	CppAD
	chol_ran_hes
	init
	pos
	ptr
	eigen
$$

$section Update Factorization Using new Matrix Values$$

$head Syntax$$
$icode%chol_ran_hes%.update(%hes_info%)%$$

$head Private$$
The $cref cholesky$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This routine updates the $cref cholesky$$ factorization
for new values in the square positive definite matrix.

$head chol_ran_hes$$
This object has prototype
$codei%
	CppAD::mixed::cholesky %chol_ran_hes%
%$$
In addition, it must have a previous call to
$cref cholesky_init$$.

$head hes_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It contains new values for the
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
we are computing the Cholesky factor of.
The $cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
must be the same as in $cref/cholesky_init/cholesky_init/hes_info/$$.
Hence, in particular, it must be in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$head ptr_$$
On input, the member variable
$codei%
	eigen_cholesky* ptr_
%$$
has been $cref/initialized/cholesky_init/$$
using the sparsity pattern for the Hessian.
Upon return, it contains the factorization
$codei%
	ptr_->factorize(%hessian%)
%$$
where $icode hessian$$ is an $code eigen_sparse$$
representation of the Hessian with values.

$end
*/
void cholesky::update(const CppAD::mixed::sparse_mat_info& hes_info)
{	assert( hes_info.row.size() == hes_info.col.size() );
	assert( hes_info.row.size() == hes_info.val.size() );
	//
	eigen_sparse hessian(n_random_, n_random_);
	for(size_t k = 0; k < hes_info.row.size(); k++)
	{	size_t r = hes_info.row[k];
		size_t c = hes_info.col[k];
		double v = hes_info.val[k];
		assert( r < n_random_ );
		assert( c < n_random_ );
		hessian.insert(r, c) = v;
	}
	// LDL^T Cholesky factorization of for specified values of the Hessian
	// f_{u,u}(theta, u)
	ptr_->factorize(hessian);
}
/*
------------------------------------------------------------------------------
$begin cholesky_logdet$$
$spell
	Cholesky
	logdet
	chol_ran_hes
	CppAD
$$

$section Compute Log Determinant for Current Cholesky Factor$$

$head Syntax$$
$icode%logdet% = %chol_ran_hes%.logdet(%n_random%)
%$$

$head Private$$
The $code cholesky$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head chol_ran_hes$$
This object has prototype
$codei%
	CppAD::mixed::cholesky %chol_ran_hes%
%$$
In addition, it must have a previous call to
$cref cholesky_update$$.

$head n_random$$
Must be the number of random effects.

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the determinant of the Hessian corresponding
to the previous call to $codei%chol_ran_hes%.factorize%$$.

$end
*/
double cholesky::logdet(size_t n_random) const
{	using Eigen::Dynamic;
    typedef Eigen::Matrix<double, Dynamic, Dynamic> dense_matrix;

	// compute the logdet( f_{u,u}(theta, u )
	dense_matrix diag = ptr_->vectorD();
	assert( diag.size() == int(n_random) );
	double logdet = 0.0;
	for(size_t j = 0; j < n_random; j++)
		logdet += log( diag(j) );

	return logdet;
}
/*
------------------------------------------------------------------------------
$begin cholesky_solve$$
$spell
	Cholesky
	logdet
	chol_ran_hes
	CppAD
	const
	eigen
$$

$section Solve Hessian Times Unknown Matrix Equals Known Matrix$$

$head Syntax$$
$icode%result% = %chol_ran_hes%.solve(%known%)
%$$

$head Private$$
The $code cholesky$$ class is an
$cref/implementation detail/cholesky/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head chol_ran_hes$$
This object has prototype
$codei%
	CppAD::mixed::cholesky %chol_ran_hes%
%$$
In addition, it must have a previous call to
$cref cholesky_update$$.

$head known$$
This argument has prototype
$codei%
	const cholesky::eigen_sparse& %known%
%$$
It is the known matrix (right hand side) in the linear equation.

$head result$$
The return value has prototype
$codei%
	cholesky::eigen_sparse& %result%
%$$
It is the solution of the equation
$codei%
	%Hessian% * %result% = %known%
%$$
where $icode Hessian$$ is the Hessian w.r.t the random effects
$latex f_{u,u} ( \theta , u )$$ corresponding to the previous call to
$cref cholesky_update$$.

$end
*/
cholesky::eigen_sparse cholesky::solve(const eigen_sparse& known) const
{	return ptr_->solve(known); }


} } // END_CPPAD_MIXED_NAMESPACE
