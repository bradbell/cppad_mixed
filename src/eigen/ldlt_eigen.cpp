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
$$

$section Update Factorization Using new Matrix Values$$

$head Syntax$$
$icode%ldlt_obj%.update(%hes_info%)%$$

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

$head pos$$
If the matrix is positive definite, $icode pos$$ is true.
Otherwise, it is false.

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
$icode%logdet% = %ldlt_obj%.logdet()%$$

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

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the determinant of the Hessian corresponding
to the previous call to $codei%ldlt_obj%.factorize%$$.

$end
*/
double ldlt_eigen::logdet(void) const
{	using Eigen::Dynamic;
    typedef Eigen::Matrix<double, Dynamic, Dynamic> dense_matrix;

	// compute the logdet( f_{u,u}(theta, u )
	dense_matrix diag = ptr_->vectorD();
	assert( diag.size() == int(n_random_) );
	double logdet = 0.0;
	for(size_t j = 0; j < n_random_; j++)
		logdet += log( diag(j) );

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
$codei%%ldlt_obj%.solve_H(%row_in%, %val_in%, %row_out%, %val_out%)%$$

$head Private$$
The $cref ldlt_cholmod$$ class is an
$cref/implementation detail/ldlt_eigen/Private/$$ and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This functions solves the linear equation
$latex H x = b$$ where $latex H$$ is the positive definite matrix
that has been factored,
$latex b$$ is a known column vector,
and $latex x$$ is unknown.

$head col_ran_hes$$
This object has prototype
$codei%
	const CppAD::mixed::ldlt_eigen %col_ran_hes%
%$$
In addition, it must have a previous call to
$cref ldlt_eigen_update$$.

$head row_in$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row_in%
%$$
It specifies the rows, in the column vector $latex b$$,
that are possibly non-zero.

$head val_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %val_in%
%$$
and it has the same size as $icode row_in$$.
It specifies the values, in the column vector $latex b$$,
that are possibly non-zero; to be specific,
for $icode%k% = 0 , %...%, %row_in%.size()-1%$$,
$codei%
	%b%[ %row_in%[%k%] ] = %val_in%[%k%]
%$$.

$head row_out$$
This argument has prototype
$codei%
	CppAD::vector<size_t>& %row_out%
%$$
On input size of $icode row_out$$ is zero.
Upon return, it contains the row indices for non-zero elements
of the solution vector $latex x$$ in increasing order; i.e.,
for $icode%k% = 1 , %...%, %row_in%.size()-1%$$,
$codei%
	%row_in%[%k-1%] < %row_out%[%k%]
%$$

$head val_out$$
This argument has prototype
$codei%
	CppAD::vector<double>& %val_in%
%$$
On input the size of $icode val_out$$ is zero.
Upon return, it has the same size at $icode row_out$$
and contains the value of the non-zero elements of the solution.
To be specific,
for $icode%k% = 0 , %...%, %row_out%.size()-1%$$,
$codei%
	%x%[ %row_out%[%k%] ] = %val_out%[%k%]
%$$

$end
*/
void ldlt_eigen::solve_H(
	const CppAD::vector<size_t>& row_in  ,
	const CppAD::vector<double>& val_in  ,
	CppAD::vector<size_t>&       row_out ,
	CppAD::vector<double>&       val_out )
{	assert( row_in.size() == val_in.size() );
	assert( row_out.size() == 0 );
	assert( val_out.size() == 0 );
	//
	eigen_sparse b(n_random_, 1);
	for(size_t k = 0; k < row_in.size(); k++)
	{	assert( row_in[k] < n_random_ );
		b.insert( row_in[k], 0 ) = val_in[k];
	}
	//
	eigen_sparse x = ptr_->solve(b);
	assert( x.outerSize() == 1 );
	assert( size_t( x.innerSize() ) == n_random_ );
	//
	typedef typename eigen_sparse::InnerIterator column_itr;
	for(column_itr itr(x, 0); itr; ++itr)
	{	assert( size_t( itr.row() ) < n_random_ );
		size_t r = size_t( itr.row() );
		double v = itr.value();
		row_out.push_back(r);
		val_out.push_back(v);
	}
}


} } // END_CPPAD_MIXED_NAMESPACE
