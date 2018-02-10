/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>

namespace {
	// --------------------------------------------------------------------
	bool check_in_limits(double lower, double x, double upper, double tol)
	{	bool flag = true;
		if( upper >= 0.0 )
			flag &= x <= (1.0 + tol) * upper;
		else
			flag &= x <= (1.0 - tol) * upper;
		//
		if( lower >= 0.0 )
			flag &= (1.0 - tol) * lower <= x;
		else
			flag &= (1.0 + tol) * lower <= x;
		//
		return flag;
	}
}

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
$begin ipopt_fixed_finalize_solution$$
$spell
	CppAD
	ran_obj
	cppad
	obj
	ipopt
	bool
	eval
	const
	obj
	ip
	cq
	namespace
	infeasibility
	doesn't
	Inf
	naninf
	std
	cout
	endl
	tol
	solution solution
$$

$section Get Solution Results$$

$head Syntax$$
$codei%finalize_solution(
	%status%, %n%, %x%, %z_L%, %z_U%, %m%, %g%,%$$
$icode%lambda%, %obj_value%, %ip_data%, %ip_cq%
)%$$

$head solution_$$
This routine checks the solution values and sets the member variable
$codei%
	CppAD::mixed fixed_solution solution_
%$$.

$head n$$
is the number of variables in the problem (dimension of x).

$head x$$
is the final value (best value found) for the primal variables
(has size $icode n$$).

$head z_L$$
is the final value for the $icode x$$ lower bound constraint multipliers
(has size $icode n$$).

$head z_U$$
is the final value for the $icode x$$ upper bound constraint multipliers
(has size $icode n$$).

$head m$$
is the number of constraints in the problem (dimension of g(x)).

$head lambda$$
is the final value for the g(x) constraint multipliers $latex \lambda$$.

$head obj_value$$
is the value of the objective f(x) at the final $icode x$$ value.

$head ip_data$$
Unspecified; i.e., not part of the Ipopt user API.

$head ip_cq$$
Unspecified; i.e., not part of the Ipopt user API.

$head status$$
These status values are in the $code Ipopt$$ namespace; e.g.,
$code SUCCESS$$ is short for $code Ipopt::SUCCESS$$:

$subhead SUCCESS$$
Algorithm terminated successfully at a locally optimal point,
satisfying the convergence tolerances (can be specified by options).

$subhead MAXITER_EXCEEDED$$
Maximum number of iterations exceeded (can be specified by an option).

$subhead CPUTIME_EXCEEDED$$
Maximum number of CPU seconds exceeded (can be specified by an option).

$subhead STOP_AT_TINY_STEP$$
Algorithm proceeds with very little progress.

$subhead STOP_AT_ACCEPTABLE_POINT$$
Algorithm stopped at a point that was converged, not to desired
tolerances, but to acceptable tolerances (see the acceptable-... options).

$subhead LOCAL_INFEASIBILITY$$
Algorithm converged to a point of local infeasibility. Problem may be
infeasible.

$subhead USER_REQUESTED_STOP$$
A user call-back function returned false, i.e.,
the user code requested a premature termination of the optimization.

$subhead DIVERGING_ITERATES$$
It seems that the iterates diverge.

$subhead RESTORATION_FAILURE$$
Restoration phase failed, algorithm doesn't know how to proceed.

$subhead ERROR_IN_STEP_COMPUTATION$$
An unrecoverable error occurred while Ipopt tried to compute
the search direction.

$subhead INVALID_NUMBER_DETECTED$$
Algorithm received an invalid number (such as NaN or Inf) from
the NLP; see also option check_derivatives_for_naninf.

$head Prototype$$
$srccode%cpp% */
void ipopt_fixed::finalize_solution(
	Ipopt::SolverReturn               status    ,  // in
	Index                             n         ,  // in
	const Number*                     x         ,  // in
	const Number*                     z_L       ,  // in
	const Number*                     z_U       ,  // in
	Index                             m         ,  // in
	const Number*                     g         ,  // in
	const Number*                     lambda    ,  // in
	Number                            obj_value ,  // in
	const Ipopt::IpoptData*           ip_data   ,  // in
	Ipopt::IpoptCalculatedQuantities* ip_cq     )  // in
/* %$$
$end
*/
{	bool ok = true;
	//
	assert( n > 0 );
	assert( size_t(n) == n_fixed_ + fix_likelihood_nabs_ );
	assert( m >= 0 );
	assert( size_t(m) == 2 * fix_likelihood_nabs_ + n_fix_con_ + n_ran_con_ );
	//
	// solution_.fixed_opt
	assert( solution_.fixed_opt.size() == 0 );
	solution_.fixed_opt.resize(n_fixed_);
	for(size_t j = 0; j < n_fixed_; j++)
		solution_.fixed_opt[j] = x[j];
	//
	// solution_.fixed_lag (see below)
	//
	// solution_.ran_con_lag
	assert( solution_.ran_con_lag.size() == 0 );
	solution_.ran_con_lag.resize(n_ran_con_);
	size_t offset = 2 * fix_likelihood_nabs_ + n_fix_con_;
	for(size_t j = 0; j < n_ran_con_; j++)
		solution_.ran_con_lag[j] = lambda[ offset + j];
	//
	// short name for fixed effects tolerance
	double tol = fixed_tolerance_;
	//
	// check that x is within its limits
	for(size_t j = 0; j < n_fixed_; j++)
	{	ok &= check_in_limits(
			fixed_lower_[j], x[j], fixed_upper_[j], 2.0 * tol
		);
	}
	//
	// check that the bound multipliers are feasible
	for(size_t j = 0; j < n_fixed_ + fix_likelihood_nabs_; j++)
	{	ok &= 0.0 <= z_L[j];
		ok &= 0.0 <= z_U[j];
	}
	//
	// fixed_opt is an alias for solution_.fixed_opt
	d_vector& fixed_opt = solution_.fixed_opt;
	//
	// fixed likelihood at the final fixed effects vector
	if( fix_likelihood_vec_tmp_.size() == 0 )
		assert( mixed_object_.fix_like_eval(fixed_opt).size() == 0 );
	else
	{	fix_likelihood_vec_tmp_ = mixed_object_.fix_like_eval(fixed_opt);
		assert( fix_likelihood_vec_tmp_.size() == 1 + fix_likelihood_nabs_ );

		// check constraints corresponding to l1 terms
		for(size_t j = 0; j < fix_likelihood_nabs_; j++)
		{	double var  = double( x[n_fixed_ + j] );
			double diff = var - fix_likelihood_vec_tmp_[j + 1];
			ok         &= scale_g_[2 * j] * diff + 1e2 * tol >= 0;
			diff        = var + fix_likelihood_vec_tmp_[j + 1];
			ok         &= scale_g_[2 * j] * diff + 1e2 * tol >= 0;
		}
	}
	//
	// explicit constraints at the final fixed effects vector
	c_vec_tmp_ = mixed_object_.fix_con_eval(fixed_opt);
	assert( c_vec_tmp_.size() == n_fix_con_ );

	// check explicit constraints and set solution_.fix_con_lag
	assert( solution_.fix_con_lag.size() == 0 );
	solution_.fix_con_lag.resize(n_fix_con_);
	offset     = 2 * fix_likelihood_nabs_;
	double inf = std::numeric_limits<double>::infinity();
	for(size_t j = 0; j < n_fix_con_; j++)
	{	// It seems from testing that Ipopt is insuring the constraint to
		// be within tol, but it seems to Brad is should be within
		// scale_g_[offset +j] * tol.
		ok &= check_in_limits(
			fix_constraint_lower_[j], c_vec_tmp_[j], fix_constraint_upper_[j],
			2.0 * tol
		);
		double lam_j  = lambda[offset + j];
		double scale = 0.0;;
		if( fix_constraint_lower_[j] != -inf )
			scale = std::fabs( fix_constraint_lower_[j] );
		if( fix_constraint_upper_[j] != inf )
			scale = std::max(scale, std::fabs( fix_constraint_upper_[j] ) );
		if( scale == 0.0 )
		{	// both limits are infinity
			lam_j = 0.0;
		}
		if( c_vec_tmp_[j] - fix_constraint_lower_[j] < tol * 10. * scale )
			lam_j = std::min(lam_j, 0.0);
		if( fix_constraint_upper_[j] - c_vec_tmp_[j] < tol * 10. * scale )
			lam_j = std::max(lam_j, 0.0);
		//
		solution_.fix_con_lag[j] = lam_j;
	}
	// Evaluate gradient of f w.r.t x
	CppAD::vector<Number> grad_f(n);
	bool new_x = true;
	eval_grad_f(n, x, new_x, grad_f.data() );

	// Evaluate gradient of g w.r.t x
	CppAD::vector<Index> iRow(nnz_jac_g_), jCol(nnz_jac_g_);
	iRow.data();
	eval_jac_g(
		n, x, new_x, m, Index(nnz_jac_g_),
		iRow.data(), jCol.data(), NULL
	);
	CppAD::vector<Number> jac_g(nnz_jac_g_);
	eval_jac_g(
		n, x, new_x, m, Index(nnz_jac_g_),
		iRow.data(), jCol.data(), jac_g.data()
	);

	// Check the partial of the Lagrangian w.r.t fixed effects
	// and set solution_.fixed_lag
	assert( solution_.fixed_lag.size() == 0 );
	solution_.fixed_lag.resize(n_fixed_);
	double average = 0.0;
	for(size_t j = 0; j < n_fixed_; j++)
	{	Number sum = grad_f[j];
		for(size_t k = 0; k < nnz_jac_g_; k++)
		{	if( jCol[k] == Index(j) )
			{	Index  i = iRow[k];
				sum   += lambda[i] * jac_g[k];
			}
		}
		// sum += z_U[j] - z_L[j]; does not work because
		// Ipopt does not seem to set z_U[j] and z_L[j] accuractely
		//
		// initialize
		solution_.fixed_lag[j] = 0.0;
		//
		// scale
		double scale = std::fabs( x[j] );
		if( fixed_lower_[j] != - inf )
			scale = std::max(scale, std::fabs( fixed_lower_[j] ) );
		if( fixed_upper_[j] != + inf )
			scale = std::max(scale, std::fabs( fixed_upper_[j] ) );
		//
		// at_lower
		bool at_lower = x[j] - fixed_lower_[j] <= scale * 10. * tol;
		// catch special case where lower and upper limits are zero
		at_lower     |= fixed_lower_[j] == fixed_upper_[j];
		at_lower     &= sum > 0.0;
		if( at_lower )
			solution_.fixed_lag[j] = - sum;
		//
		// at_upper
		bool at_upper = fixed_upper_[j] - x[j] <= scale * 10. * tol;
		at_upper     |= fixed_lower_[j] == fixed_upper_[j];
		at_upper     &= sum < 0.0;
		if( at_upper )
			solution_.fixed_lag[j] = - sum;
		//
		if( ! (at_lower || at_upper) )
		{	double check;
			check = std::fabs(sum);
			if( sum >= 0.0  && fixed_lower_[j] > -inf )
				check = sum * (x[j] - fixed_lower_[j]);
			else if( sum <= 0.0 && fixed_upper_[j] < +inf )
				check = - sum * (fixed_upper_[j] - x[j]);
			assert( check >= 0.0 );
			//
			average += check / double(n_fixed_);
		}
	}
	ok &= average * scale_f_ <= 10. * tol;

	// Check the partial of the Lagrangian w.r.t auxillary variables
	average = 0.0;
	for(size_t j = n_fixed_; j < n_fixed_ + fix_likelihood_nabs_; j++)
	{	Number sum = grad_f[j];
		for(size_t k = 0; k < nnz_jac_g_; k++)
		{	if( jCol[k] == Index(j) )
			{	Index  i = iRow[k];
				sum   += lambda[i] * jac_g[k];
			}
		}
		average += std::fabs(sum) / double(n_fixed_);
	}
	ok &= average * scale_f_ <= tol;

	// set member variable finalize_solution_ok_
	finalize_solution_ok_ = ok;
}
} } // END_CPPAD_MIXED_NAMESPACE
