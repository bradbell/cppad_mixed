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
cholesky::cholesky(void)
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
$icode%chol_ran_hes%.init(%n_fixed%, %n_random%, %row%, %col%)
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

$head n_fixed$$
This argument has prototype
$codei%
	size_t %n_fixed%
%$$
and is then number of fixed effects.

$head n_random$$
This argument has prototype
$codei%
	size_t %n_random%
%$$
and is then number of random effects.

$head row$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %row%
%$$
These are the possibly non-zero row indices in the sparsity pattern
of the lower triangle of the Hessian w.r.t just the random effects.
For $icode%k% = 0 , %...% , %row%.size()-1%$$,
$codei%
	%n_fixed% <= %row%[%k%] < %n_fixed% + %n_random%
%$$
The reason for the offset is that these indices are relative to both
the fixed and random effects and the fixed effects come before the
random effects.

$head col$$
This argument has prototype
$codei%
	const CppAD::vector<size_t>& %col%
%$$
These are the possibly non-zero column indices in the sparsity pattern
of the lower triangle of the Hessian.
It must have the same size as $icode row$$ and
for $icode%k% = 0 , %...% , %col%.size()-1%$$,
$codei%
	%n_fixed% <= %col%[%k%] <= %row%[%k%]
%$$

$end
*/
void cholesky::init(
	size_t                       n_fixed  ,
	size_t                       n_random ,
	const CppAD::vector<size_t>& row      ,
	const CppAD::vector<size_t>& col      )
{	double not_used = 1.0;

	Eigen::SparseMatrix<double> hessian_pattern(n_random, n_random);
	assert( row.size() == col.size() );
	for(size_t k = 0; k < row.size(); k++)
	{
		assert( n_fixed <= row[k] && row[k] < n_fixed + n_random );
		assert( n_fixed <= col[k] && col[k] <= row[k] );
		hessian_pattern.insert(row[k] - n_fixed, col[k] - n_fixed) = not_used;
	}
	// analyze the pattern for an LDL^T Cholesky factorization of
	// f_{u,u}(theta, u)
	ptr_->analyzePattern(hessian_pattern);
}
/*
-----------------------------------------------------------------------------
$begin cholesky_factorize$$
$spell
	Cholesky
	chol_ran_hes
	CppAD
	const
	Taylor
	init
$$

$section Compute Cholesky Factor for Specific Fixed and Random Effects$$

$head Syntax$$
$icode%chol_ran_hes%.factorize(
	%n_fixed%, %n_random%, %row%, %col%, %both%, %ran_hes_fun%
)
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
$cref cholesky_init$$.

$head n_fixed$$
Must be the same as $icode n_fixed$$ in previous call to
$cref/cholesky_init/cholesky_init/n_fixed/$$.

$head n_random$$
Must be the same as $icode n_random$$ in previous call to
$cref/cholesky_init/cholesky_init/n_random/$$.

$head row$$
Must be the same as $icode row$$ in previous call to
$cref/cholesky_init/cholesky_init/row/$$.

$head col$$
Must be the same as $icode col$$ in previous call to
$cref/cholesky_init/cholesky_init/col/$$.

$head both$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %both%
%$$
and has size $icode%n_fixed% + %n_random%$$.
This is the values of the fixed and random effects at which the Hessian
is being computed and factored.
The fixed effects come first and then the random effects.

$head ran_hes_fun$$
This argument has prototype
$codei%
	CppAD::ADFun<double>& %ran_hes_fun%
%$$
It has $icode%ran_hes_fun%.Domain() = %n_fixed% + %n_random%$$
and $icode%ran_hes_fun%.Range() = %row%.size()%$$
The function call
$codei%
	%val% = %ran_hes_fun%.Forward(0, %both%)
%$$
is used to compute the values of the Hessian.
Thus, upon return, the first order Taylor coefficient for the corresponding
fixed and random effects are stored in $icode ran_hes_fun$$.

$end
*/

void cholesky::factorize(
	size_t                       n_fixed      ,
	size_t                       n_random     ,
	const CppAD::vector<size_t>& row          ,
	const CppAD::vector<size_t>& col          ,
	const CppAD::vector<double>& both         ,
	CppAD::ADFun<double>&        ran_hes_fun  )
{
	size_t K = row.size();
	CppAD::vector<double> val(K);
	val = ran_hes_fun.Forward(0, both);

	Eigen::SparseMatrix<double> hessian_value(n_random, n_random);
	assert( row.size() == col.size() );
	for(size_t k = 0; k < row.size(); k++)
	{	assert( n_fixed <= row[k] && row[k] < n_fixed + n_random );
		assert( n_fixed <= col[k] && col[k] <= row[k] );
		hessian_value.insert(row[k] - n_fixed, col[k] - n_fixed) = val[k];
	}
	// LDL^T Cholesky factorization of for specified values of the Hessian
	// f_{u,u}(theta, u)
	ptr_->factorize(hessian_value);
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
$cref cholesky_factorize$$.

$head n_random$$
Must be the same as $icode n_random$$ in previous call to
$cref/cholesky_factorize/cholesky_factorize/n_random/$$.

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
$cref cholesky_factorize$$.

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
$cref cholesky_factorize$$.

$end
*/
cholesky::eigen_sparse cholesky::solve(const eigen_sparse& known) const
{	return ptr_->solve(known); }


} } // END_CPPAD_MIXED_NAMESPACE
