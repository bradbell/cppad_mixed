// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: Estimating Disease Rates as Functions of Age and Time
          Copyright (C) 2014-15 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/chol_hes_ran.hpp>
# include <Eigen/Sparse>

/*
$begin chol_hes_ran$$
$spell
	CppAD
	init_chol_hes
	CppAD
	cppad
	const
	Cholesky
	namespace
	Eigen
	hpp
	Simplicial
	triangular
	dismod
	logdet
$$

$section Sparse Cholesky Factorization of Hessian w.r.t Random Effects$$

$head Syntax$$
$codei%CppAD::mixed::analyze_chol_hes_ran(%n_fixed%, %n_random%, %row%, %col%)
%$$
$codei%CppAD::mixed::factorize_chol_hes_ran(
	%n_fixed%, %n_random%, %row%, %col%, %both%, %hessian%
)%$$

$head Purpose$$
$code chol_hes_ran_$$ should be a $cref private$$ member variable,
but it is instead a static in the $code cppad_mixed$$ namespace
so that the warnings that Eigen generates
do not need to be suppressed by all the routines that include
$code cppad_mixed/cppad_mixed.hpp$$.

$head chol_hes_ran_$$
This variable has prototype
$codep */
	namespace CppAD { namespace mixed {
		Eigen::SimplicialLDLT<
			Eigen::SparseMatrix<double> , Eigen::Lower
		> chol_hes_ran_;
	} }
/* $$
This is lower triangular Cholesky factorization of the Hessian of the
$cref/random likelihood
	/theory
	/Random Likelihood, f(theta, u)
/$$
with respect to the random effects; i.e.
$latex f_{uu}^{(2)} ( \theta , u )$$.


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
For $icode%k% = 0 , %...% , %row%.size()%$$,
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
for $icode%k% = 0 , %...% , %col%.size()%$$,
$codei%
	%n_fixed% <= %col%[%k%] <= %row%[%k%]
%$$

$head both$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %both%
%$$
and has size $icode%n_fixed% + %n_random%$$.
This is the values of the fixed and random effects at which the Hessian
is being computed and factored.

$head hessian$$
This argument has prototype
$codei%
	CppAD::ADFun<double>& %hessian%
%$$
It has $icode%hessian%.Domain() = %n_fixed% + %n_random%$$
and $icode%hessian%.Range() = %row%.size()%$$
The function call
$codei%
	%val% = %hessian%.Forward(0, %both%)
%$$
is used to compute the values of the Hessian.
To be specific,
for $icode%k% = 0 , %...% , %col%.size()%%$$,
$icode%val%[%k%]%$$ is the value of the Hessian
at row index $icode%row%[%k%]%$$
and column index $icode%col%[%k%]%$$.

$head analyze_chol_hes_ran$$
The input value of this factorization does not matter.
Upon return, the sparsity pattern has been analyzed; i.e.,
$codei%
	chol_hes_ran_.analyzePattern(%hessian_pattern%)
%$$
has been called with the pattern corresponding to $icode row$$ and
$icode col$$.

$head factorize_chol_hes_ran$$
The $icode row$$ and $icode col$$ values must be the same as for
the previous call to
$cref/analyze_chol_hes_ran/chol_hes_ran/analyze_chol_hes_ran/$$.
Upon return, the Hessian has been evaluated and then factorized; i.e.,
$codei%
	chol_hes_ran_.factorize(%hessian_value%)
%$$
has been called with the values corresponding to the Hessian at the
fixed and random effects specified by $icode both$$.

$head logdet$$
This return value has prototype
$codei%
	double %logdet%
%$$
Is the log of the determinant of the Hessian corresponding
to the previous call to
$cref/factorize_chol_hes_ran/chol_hes_ran/factorize_chol_hes_ran/$$.

$end
*/

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

void analyze_chol_hes_ran(
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
	// f_{uu}^{(2)}(theta, u)
	chol_hes_ran_.analyzePattern(hessian_pattern);
}

void factorize_chol_hes_ran(
	size_t                       n_fixed  ,
	size_t                       n_random ,
	const CppAD::vector<size_t>& row      ,
	const CppAD::vector<size_t>& col      ,
	const CppAD::vector<double>& both     ,
	CppAD::ADFun<double>&        hessian  )
{
	size_t K = row.size();
	CppAD::vector<double> val(K);
	val = hessian.Forward(0, both);

	Eigen::SparseMatrix<double> hessian_value(n_random, n_random);
	assert( row.size() == col.size() );
	for(size_t k = 0; k < row.size(); k++)
	{	assert( n_fixed <= row[k] && row[k] < n_fixed + n_random );
		assert( n_fixed <= col[k] && col[k] <= row[k] );
		hessian_value.insert(row[k] - n_fixed, col[k] - n_fixed) = val[k];
	}
	// LDL^T Cholesky factorization of for specified values of the Hessian
	// f_{uu}^{(2)}(theta, u)
	chol_hes_ran_.factorize(hessian_value);
}

double logdet_chol_hes_ran(size_t n_random)
{	using Eigen::Dynamic;
    typedef Eigen::Matrix<double, Dynamic, Dynamic> dense_matrix;

	// compute the logdet( f_{uu}^{(2)}(theta, u )
	dense_matrix diag = chol_hes_ran_.vectorD();
	assert( diag.size() == int(n_random) );
	double logdet = 0.0;
	for(size_t j = 0; j < n_random; j++)
		logdet += log( diag(j) );

	return logdet;
}


} } // END_CPPAD_MIXED_NAMESPACE
