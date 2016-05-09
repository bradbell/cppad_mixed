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
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/triple2eigen.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ldlt_eigen_ctor$$
$spell
	ldlt_obj
	eigen
	ptr
	CppAD
$$

$section Eigen LDLT Constructor$$

$head Syntax$$
$codei%CppAD::mixed::ldlt_eigen %ldlt_obj%(%n_random%)%$$


$head Private$$
The $cref ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head n_random$$
The argument $icode n_random$$ has prototype
$codei%
	size_t %n_random%
%$$
It is the number of rows in the symmetric matrix we will compute factor.
The member variable $code n_random_$$ is set to this value.

$head ptr_$$
This member variable points to a newly constructed
$code eigen_ldlt$$ object.

$end
*/

ldlt_eigen::ldlt_eigen(size_t n_random)
: n_random_(n_random)
{	ptr_ = new eigen_ldlt; }

// destructor
ldlt_eigen::~ldlt_eigen(void)
{	delete ptr_; }

/*
------------------------------------------------------------------------------
$begin ldlt_eigen_init$$
$spell
	pos
	ldlt_eigen
	ldlt_obj
	CppAD
	const
	init
	hes
$$

$section Initialize LDLT Factor for a Specific Sparsity Pattern$$

$head Syntax$$
$icode%pos% = %ldlt_obj%.init(%hes_info%)
%$$

$head Private$$
The $code ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head ldlt_obj$$
This object has prototype
$codei%
	CppAD::mixed::ldlt_eigen %ldlt_obj%
%$$

$head hes_info$$
This argument had prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It is a
$cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$ for the
matrix that we will compute the LDLT factor of.
It is in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$end
*/
void ldlt_eigen::init( const CppAD::mixed::sparse_mat_info& hes_info )
{	assert( hes_info.row.size() == hes_info.col.size() );
	CppAD::vector<double> not_used(0);
	//
	eigen_sparse hessian_pattern = CppAD::mixed::triple2eigen(
			n_random_           ,
			n_random_           ,
			hes_info.row        ,
			hes_info.col        ,
			not_used
	);
	// analyze the pattern for an LDLT factorization of
	// f_{u,u}(theta, u)
	ptr_->analyzePattern(hessian_pattern);
}
/*
------------------------------------------------------------------------------
$begin ldlt_eigen_update$$
$spell
	ldlt_eigen
	xam
	const
	CppAD
	ldlt_obj
	init
	pos
	ptr
	eigen
	hes
	bool
$$

$section Update Factorization Using new Matrix Values$$

$head Syntax$$
$icode%flag% = ldlt_obj%.update(%hes_info%)%$$

$head Private$$
The $cref ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This routine updates the $cref ldlt_eigen$$ factorization
for new values in the square positive definite matrix.

$head ldlt_obj$$
This object has prototype
$codei%
	CppAD::mixed::ldlt_eigen %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_eigen_init$$.

$head hes_info$$
This argument has prototype
$codei%
	const CppAD::mixed::sparse_mat_info& %hes_info%
%$$
It contains new values for the
$cref/sparse matrix/sparse_mat_info/Notation/Sparse Matrix/$$
we are computing the LDLT factor of.
The $cref/sparsity pattern/sparse_mat_info/Notation/Sparsity Pattern/$$
must be the same as in $cref/ldlt_eigen_init/ldlt_eigen_init/hes_info/$$.
Hence, in particular, it must be in
$cref/column major/sparse_mat_info/Notation/Column Major Order/$$ order
and
$cref/lower triangular/sparse_mat_info/Notation/Lower Triangular/$$.

$head ptr_$$
On input, the member variable
$codei%
	eigen_ldlt* ptr_
%$$
has been $cref/initialized/ldlt_eigen_init/$$
using the sparsity pattern for the Hessian.
Upon return, it contains the factorization
$codei%
	ptr_->factorize(%hessian%)
%$$
where $icode hessian$$ is an $code eigen_sparse$$
representation of the Hessian with values.

$head flag$$
The return value has prototype
$codei%
	bool %flag%
%$$
It is true for success and false for numerical issues.

$end
*/
bool ldlt_eigen::update(const CppAD::mixed::sparse_mat_info& hes_info)
{	assert( hes_info.row.size() == hes_info.col.size() );
	assert( hes_info.row.size() == hes_info.val.size() );
	//
	eigen_sparse hessian = CppAD::mixed::triple2eigen(
		n_random_      ,
		n_random_      ,
		hes_info.row   ,
		hes_info.col   ,
		hes_info.val
	);
	// LDLT factorization of for specified values of the Hessian
	// f_{u,u}(theta, u)
	ptr_->factorize(hessian);
	//
	if( ptr_->info() == Eigen::NumericalIssue )
		return false;
	assert( ptr_->info() == Eigen::Success );
	return true;
}
/*
------------------------------------------------------------------------------
$begin ldlt_eigen_logdet$$
$spell
	ldlt_eigen
	logdet
	ldlt_obj
	CppAD
	const
$$

$section Compute Log Determinant for Current LDLT Factor$$

$head Syntax$$
$icode%logdet% = %ldlt_obj%.logdet(%sign%)%$$

$head Private$$
The $code ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_eigen %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_eigen_update$$.

$head sign$$
This argument has prototype
$codei%
	int& %sign%
%$$
Its input value does no matter,
upon return it is $code +1$$ if the determinant is positive,
$code 0$$ if the determinant is zero,
and $code -1$$ if the determinant is negative.

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the absolute value of the determinant corresponding
to the previous call to $codei%ldlt_obj%.factorize%$$.

$end
*/
double ldlt_eigen::logdet(int& sign) const
{	using Eigen::Dynamic;
    typedef Eigen::Matrix<double, Dynamic, Dynamic> dense_matrix;

	// compute the logdet( f_{u,u}(theta, u )
	dense_matrix diag = ptr_->vectorD();
	assert( diag.size() == int(n_random_) );
	double logdet = 0.0;
	sign = 1;
	for(size_t j = 0; j < n_random_; j++)
	{	if( diag(j) == 0.0 )
		{	sign = 0;
			return - std::numeric_limits<double>::infinity();
		}
		if( diag(j) < 0.0 )
			sign = - sign;
		logdet += log( std::fabs( diag(j) ) );
	}

	return logdet;
}
/*
-----------------------------------------------------------------------------
$begin ldlt_eigen_solve_H$$
$spell
	ldlt_eigen
	ldlt_obj
	cholesky
	const
	xam
	CppAD
	hes
$$

$section Solve Linear Equations Using LDLT Factor$$

$head Syntax$$
$codei%%ldlt_obj%.solve_H(%row%, %val_in%, %val_out%)%$$

$head Private$$
The $cref ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function solves the linear equation
$latex H x = b$$ where $latex H$$ is the positive definite matrix
that has been factored,
$latex b$$ is a known column vector,
and $latex x$$ is unknown.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_eigen %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_eigen_update$$.

$head row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
It contains all of the rows of column vector $latex b$$ that are
non-zero and the rows of the column vector $icode x$$
that are desired.
These values are in strictly increasing order; i.e.,
$codei%
	%row%[%k%] < %row%[%k%+1]
%$$
It follows that $icode%row%.size()%$$ is less than or equal
$cref/n_random/ldlt_eigen_ctor/n_random/$$.

$head val_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_in%
%$$
and it has the same size as $icode row$$.
It specifies the values in the column vector $latex b$$
for each of the corresponding rows; i.e.,
for $icode%k% = 0 , %...%, %row%.size()-1%$$,
$codei%
	%b%[ %row%[%k%] ] = %val_in%[%k%]
%$$.

$head val_out$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_out%
%$$
and it has the same size as $icode row$$.
On input, the value of its elements do not matter.
Upon return, it contains the values in the column vector $latex b$$
for each of the corresponding rows; i.e.,
for $icode%k% = 0 , %...%, %row%.size()-1%$$,
$codei%
	%x%[ %row%[%k%] ] = %val_out%[%k%]
%$$.


$end
*/
void ldlt_eigen::solve_H(
	const CppAD::vector<size_t>& row     ,
	const CppAD::vector<double>& val_in  ,
	CppAD::vector<double>&       val_out )
{	assert( row.size() == val_in.size() );
	assert( row.size() == val_out.size() );
	//
	eigen_sparse b(n_random_, 1);
	for(size_t k = 0; k < row.size(); k++)
	{	assert( row[k] < n_random_ );
		b.insert( row[k], 0 ) = val_in[k];
		val_out[k] = 0.0;
	}
	//
	eigen_sparse x = ptr_->solve(b);
	assert( x.outerSize() == 1 );
	assert( size_t( x.innerSize() ) == n_random_ );
	//
# ifndef NDEBUG
	int previous_row = -1;
# endif
	size_t k = 0;
	typedef typename eigen_sparse::InnerIterator column_itr;
	for(column_itr itr(x, 0); itr; ++itr)
	{	assert( size_t( itr.row() ) < n_random_ );
		assert( previous_row < itr.row() );
# ifndef NDEBUG
		previous_row = itr.row();
# endif
		size_t r = size_t( itr.row() );
		while( k < row.size() && row[k] < r )
			k++;
		if( k < row.size() && row[k] == r )
			val_out[k] = itr.value();
	}
}
/*
$begin ldlt_eigen_sim_cov$$
$spell
	ldlt_obj
	sim_cov
	const
	eigen
	bool
	xam
	CppAD
$$

$section Simulations with Covariance Corresponding to Factored Matrix$$

$head Syntax$$
$icode%ok% = %ldlt_obj%.sim_cov(%w%, %v%)%$$

$head Private$$
The $cref ldlt_eigen$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This function simulates a normal random vector with mean zero
and covariance $latex H^{-1}$$
where $latex H$$ is the symmetric matrix that has been factored as
$latex \[
	L D L^\R{T} = P H P^\R{T}
\] $$
where $latex L$$ is lower triangular, $latex D$$ is diagonal,
and $latex P$$ is a permutation matrix.

$head ldlt_obj$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_eigen %ldlt_obj%
%$$
In addition, it must have a previous call to
$cref ldlt_eigen_update$$.

$head w$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %w%
%$$
and its size is equal to the number of rows in $latex H$$.

$head v$$
This argument has prototype
$codei%
	CppAD::vector<double>& %v%
%$$
and its size is equal to the number of rows in $latex H$$.
The input value of its elements does not matter.
Upon return
$latex \[
	v = P^\R{T} L^{-\R{T}} D^{-1/2} w
\] $$
If $latex w$$ is mean zero, variance identity white noise,
$latex w \sim \B{N} ( 0 , I )$$,
then $latex v$$ will be mean zero and variance $latex H^{-1}$$,
$latex v \sim \B{N} ( 0 , H^{-1} )$$; see
$cref/sparse observed information/theory/Sparse Observed Information/$$.

$head ok$$
The return value has prototype
$codei%
	bool %ok%
%$$
It is true if all the elements of $latex D$$ are greater then zero; i.e.,
if $latex H$$ is positive definite.
If $icode ok$$ is false,
$latex H$$ is not positive definite and the values in $latex v$$ are
the same as their input values.

$head Example$$
The file $cref/ldlt_eigen_xam.cpp/ldlt_eigen_xam.cpp/sim_cov/$$ contains an
example and test that uses this function.

$end
*/
bool ldlt_eigen::sim_cov(
	const CppAD::vector<double>& w  ,
	CppAD::vector<double>&       v  )
{	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> column_vector;
	//
	// set b = w
	column_vector b(n_random_);
	for(size_t i = 0; i < n_random_; i++)
		b[i] = w[i];
	//
	// set b = D^{-1/2} w
	column_vector diag = ptr_->vectorD();
	for(size_t i = 0; i < n_random_; i++)
	{	if( diag[i] <= 0.0 )
			return false;
		b[i] = b[i] / std::sqrt( diag[i] );
	}
	//
	// set b = L^{-T} * D^{-1/2} w
	b = ptr_->matrixU().solve(b);
	//
	// set b = P^T L^{-T} * D^{-1/2} w
	b = ptr_->permutationP().transpose() * b;
	//
	// return v
	for(size_t i = 0 ; i < n_random_; i++)
		v[i] = b[i];
	//
	return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
