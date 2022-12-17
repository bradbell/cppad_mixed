/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/exception.hpp>
/*
$begin logdet_jac$$
$spell
	ldlt
	Jacobian
	jac
	CppAD
	cppad
	hes
	vec
	const
	Cpp
	xam
	logdet
	cholesky
$$

$nospell
$bold This is old cppad_mixed documentation:$$ Here is a link to its
$href%http://bradbell.github.io/cppad_mixed%current documentation%$$.
$$
$section Jacobian of Log Determinant of Hessian w.r.t. Random Effects$$

$head Syntax$$
$icode%mixed_object%.logdet_jac(
	%fixed_vec%, %random_vec%, %logdet_fix%, %logdet_ran%
)%$$

$head Private$$
This $code cppad_mixed$$ is a $cref private_base_class$$ member function.

$head Purpose$$
This routine computes the Jacobian of the log determinant
of the Hessian of the random likelihood
$cref/f(theta, u)
	/theory/
	Random Likelihood, f(theta, u)
/$$
with respect to the random effects vector $latex u$$.
To be specific, it computes both
$latex \[
	\partial_\theta \log \det [ f_{u,u} ( \theta, u ) ]
	\; \R{and} \;
	\partial_u \log \det [ f_{u,u} ( \theta, u ) ]
\] $$

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head ldlt_ran_hes_$$
It is assumed that the member variable
$codei%
	CPPAD_MIXED_LDLT ldlt_ran_hes_
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

$head logdet_fix$$
This argument has prototype
$codei%
	CppAD::vector<double>& %logdet_fix%
%$$
Its input size must be equal to $code n_fixed_$$.
Upon return, it contains the value of the derivative w.r.t
the fixed effects.
$codei%ran_hes_.col[%k%] - n_fixed_%$$ of the Hessian.

$head logdet_ran$$
This argument has prototype
$codei%
	CppAD::vector<double>& %logdet_ran%
%$$
Its input size must be equal to $code n_random_$$.
Upon return, it contains the value of the derivative w.r.t
the random effects.

$children%
	example/private/logdet_jac.cpp
%$$
$head Example$$
The file $cref logdet_jac.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/
// ----------------------------------------------------------------------------
void cppad_mixed::logdet_jac(
	const d_vector& fixed_vec  ,
	const d_vector& random_vec ,
	d_vector&       logdet_fix ,
	d_vector&       logdet_ran )
{	assert( init_ran_hes_done_ );
	//
	assert( fixed_vec.size() == n_fixed_ );
	assert( random_vec.size() == n_random_ );
	assert( logdet_fix.size() == n_fixed_ );
	assert( logdet_ran.size() == n_random_ );
	//
	// compute the inverse where Hessian is possibly non-zero
	CppAD::mixed::sparse_mat_info weight_info;
	size_t K = ran_hes_uu_rcv_.nnz();
	weight_info.row.resize(K);
	weight_info.col.resize(K);
	weight_info.val.resize(K);
	for(size_t k = 0; k < K; k++)
	{	size_t r = ran_hes_uu_rcv_.row()[k];
		size_t c = ran_hes_uu_rcv_.col()[k];
		assert(r < n_random_);
		assert(c < n_random_);
		weight_info.row[k] = r;
		weight_info.col[k] = c;
	}
	ldlt_ran_hes_.inv(
			weight_info.row,
			weight_info.col,
			weight_info.val
	);
	// must weight the off diagonal elements twice
	// (to account for upper diagonal entry which is not present)
	for(size_t k = 0; k < K; k++)
	{	if( weight_info.row[k] != weight_info.col[k] )
			weight_info.val[k] *= 2.0;
	}
	d_vector dw(n_fixed_ + n_random_);
	dw = ran_hes_fun_.Reverse(1, weight_info.val);
	if( CppAD::hasnan( dw ) ) throw CppAD::mixed::exception(
		"logdet_jac", "result has a nan"
	);
	//
	// split out fixed and random parts of the derivative
	unpack(logdet_fix, logdet_ran, dw);
	//
	return;
}
