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
$begin optimize_fixed$$
$spell
	CppAD
	cppad
	ipopt
	xam
	vec
	const
	CppAD
	std
	inf
	iter
$$

$section Optimize Fixed Effects$$

$head Syntax$$
$icode%solution% =%$$
$icode%mixed_object%.optimize_fixed(
	%fixed_options%,
	%random_options%,
	%fixed_lower%,
	%fixed_upper%,
	%fix_constraint_lower%,
	%fix_constraint_upper%,
	%fixed_in%,
	%random_lower%,
	%random_upper%,
	%random_in%
)%$$

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head Purpose$$
This routine maximizes the total objective
$cref/L(theta)/theory/Objective/Total Objective, L(theta)/$$.

$head inf$$
The value $code inf$$ below refers to
$codei%
	std::numeric_limits<double>::infinity()
%$$

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head fixed_options$$
This argument has prototype
$codei%
	const std::string& %fixed_options%
%$$
and is the $cref ipopt_options$$ for optimizing the fixed effects.
If $cref/quasi_fixed/derived_ctor/quasi_fixed/$$
is true,
the following changes are made to the standard Ipopt options specification:
$list number$$
If specified, the $icode derivative_test$$ must be $code none$$ or
$code first-order$$.
$lnext
The value used for $icode hessian_approximation$$ is $code limit-memory$$
and cannot be specified differently in $icode fixed_options$$.
$lnext
The default value used for $icode limited_memory_max_history$$ is 30.
$lnext
If $icode max_iter$$ is zero,
$icode%solution%.fixed_opt = %fixed_in%$$
(this is not necessarily true for Ipopt).
$lend


$head random_options$$
This argument has prototype
$codei%
	const std::string& %random_options%
%$$
and is the $cref ipopt_options$$ for optimizing the random effects.


$head fixed_lower$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_lower%
%$$
It has size $cref/n_fixed/derived_ctor/n_fixed/$$ and
specifies the lower limits for the
$fixed_effects/cppad_mixed/Fixed Effects, theta/$$.
Note that minus infinity is used for no lower limit.

$head fixed_upper$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_upper%
%$$
It has size $icode n_fixed$$ and specifies the upper limits for the
$fixed_effects/cppad_mixed/Fixed Effects, theta/$$.
Note that plus infinity is used for no upper limit.

$head fix_constraint_lower$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fix_constraint_lower%
%$$
it has size $icode n_fixed$$ and specifies the lower limits for the
$cref/fixed constraints/fix_constraint/$$.
Note that minus infinity is used for no lower limit.

$head fix_constraint_upper$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fix_constraint_upper%
%$$
it specifies the upper limits for the
$cref/fixed constraints/fix_constraint/$$.
Note that plus infinity is used for no upper limit.

$head fixed_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_in%
%$$
It specifies the initial value for the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$ during the optimization process.
It must hold for each $icode j$$ that
$codei%
	%fixed_lower%[%j%] <= %fixed_in%[%j%] <= %fixed_upper%[%j%]
%$$


$head random_lower$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_lower%
%$$
It must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the lower limits for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.
This may be useful to keep the random effects
out of regions of numerical instability.
On the other hand, the calculation of the
$cref/derivative of optimal random effects
	/theory
	/Derivative of Optimal Random Effects
/$$
$latex \hat{u}_\theta ( \theta )$$ will not be correct when these constraints
are active (and this could have adverse effects on the optimization).
The value minus infinity can be used to specify no lower limit.

$head random_upper$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_upper%
%$$
It must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the upper limits for the optimization of the random effect.
The value plus infinity can be used to specify no lower limit.

$head random_in$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %random_in%
%$$
It must have size equal to
$cref/n_random/derived_ctor/n_random/$$ and
specifies the initial value used for the optimization of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$ vector $latex u$$.
It must hold that
$codei%
	%random_lower%[%i%] <= %random_in%[%i%] <= %random_upper%[%i%]
%$$
for each valid index $icode i$$.

$head solution$$
The return value has prototype
$codei%
	CppAD::mixed::fixed_solution %solution%
%$$
It is the solution (obtained by optimization) of the
fixed effects vector and its Lagrange multipliers; see
$cref fixed_solution$$.

$head Laplace Approximation$$
The $cref theory$$ for the
Laplace approximation optimization only includes the case where
the $cref/random likelihood/ran_likelihood/$$ is smooth.

$comment ipopt_options is also used by optimize_random$$
$children%example/user/optimize_fixed_xam.cpp
	%src/ipopt_options.omh
%$$
$head Example$$
The file $cref optimize_fixed_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$head ipopt_fixed$$
The  class $cref ipopt_fixed$$ is used by $code optimize_fixed$$
to optimize the fixed effects.
It's specifications are not part of the $cref cppad_mixed$$ public interface.

$end
------------------------------------------------------------------------------
*/
# include <coin/IpIpoptApplication.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/ipopt_fixed.hpp>



CppAD::mixed::fixed_solution cppad_mixed::optimize_fixed(
	const std::string& fixed_options     ,
	const std::string& random_options    ,
	const d_vector&    fixed_lower       ,
	const d_vector&    fixed_upper       ,
	const d_vector&    fix_constraint_lower  ,
	const d_vector&    fix_constraint_upper  ,
	const d_vector&    fixed_in          ,
	const d_vector&    random_lower      ,
	const d_vector&    random_upper      ,
	const d_vector&    random_in         )
{	bool ok = true;
	//
	// fixed_(lower, upper, in)
	assert( fixed_lower.size() == n_fixed_ );
	assert( fixed_upper.size() == n_fixed_ );
	assert( fixed_in.size()    == n_fixed_ );
	// fix_constraint(lower and upper)
	assert( fix_constraint_lower.size() == fix_con_fun_.Range() );
	assert( fix_constraint_upper.size() == fix_con_fun_.Range() );
	// random_(lower, upper, in)
	assert( random_lower.size() == n_random_ );
	assert( random_upper.size() == n_random_ );
	assert( random_in.size() == n_random_ );
	//
	using Ipopt::SmartPtr;
	// make sure initialize has been called
	if( ! initialize_done_ )
	{	std::string error_message =
		"cppad_mixed::initialize was not called before optimize_fixed";
		fatal_error(error_message);
	}

# ifndef NDEBUG
	for(size_t j = 0; j < n_fixed_; j++)
	{	assert( fixed_lower[j] <= fixed_in[j] );
		assert( fixed_in[j]    <= fixed_upper[j] );
	}
# endif

	// create a reference to this object
	cppad_mixed& mixed_object(*this);

	// Create an instance of an IpoptApplication
	SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

	if( quasi_fixed_ )
	{	// special defaults settings
		app->Options()->SetStringValue(
			"hessian_approximation", "limited-memory");
		app->Options()->SetIntegerValue(
			"limited_memory_max_history", 30);
	}
	// Set options for optimization of the fixed effects
	const std::string& options = fixed_options;
	size_t begin_1, end_1, begin_2, end_2, begin_3, end_3;
	begin_1     = 0;
	while( options[begin_1] == ' ')
		begin_1++;
	while( begin_1 < options.size() )
	{	// split this line into tokens
		end_1   = options.find_first_of(" \n", begin_1);
		begin_2 = end_1;
		while( options[begin_2] == ' ')
			begin_2++;
		end_2   = options.find_first_of(" \n", begin_2);
		begin_3 = end_2;
		while( options[begin_3] == ' ')
			begin_3++;
		end_3   = options.find_first_of(" \n", begin_3);

		// check for three non-empty tokens
		assert( end_3 != std::string::npos );
		assert( begin_1 < end_1 && end_1 <= begin_2 );
		assert( begin_2 < end_2 && end_2 <= begin_3 );
		assert( begin_3 < end_3 );

		// get the three tokens
		std::string tok_1 = options.substr(begin_1, end_1 - begin_1);
		std::string tok_2 = options.substr(begin_2, end_2 - begin_2);
		std::string tok_3 = options.substr(begin_3, end_3 - begin_3);
		// switch on option type
		if ( tok_1 == "String" )
		{	app->Options()->SetStringValue(tok_2.c_str(), tok_3.c_str());
			if( quasi_fixed_ )
			{	ok = true;
				if( tok_2 == "hessian_approximation" )
					ok &= tok_3 == "limited-memory";
				if( tok_2 == "derivative_test" )
					ok &= tok_3 == "none" ||tok_3 == "first-order";
				if( ! ok )
				{	std::string msg = "cppad_mixed: constructed with";
					msg += " quasi_fixed true so cannot have ";
					msg += tok_2 + " equal to " + tok_3;
					fatal_error(msg);
				}
			}
		}
		else if ( tok_1 == "Numeric" )
		{	Ipopt::Number value = std::atof( tok_3.c_str() );
			app->Options()->SetNumericValue(tok_2.c_str(), value);
		}
		else if ( tok_1 == "Integer" )
		{	Ipopt::Index value = std::atoi( tok_3.c_str() );
			app->Options()->SetIntegerValue(tok_2.c_str(), value);
		}
		else assert(false);
		//
		// next line
		begin_1 = end_3;
		while( options[begin_1] == ' ' || options[begin_1] == '\n' )
			begin_1++;
	}
	// get the tolerance settting for the fixed effects optimization
	std::string tag    = "tol";
	std::string prefix = "";
	double fixed_tolerance;
	app->Options()->GetNumericValue(tag, fixed_tolerance, prefix);
	// get the maximum number of iterations
	tag    = "max_iter";
	prefix = "";
	Ipopt::Index fixed_max_iter;
	app->Options()->GetIntegerValue(tag, fixed_max_iter, prefix);
	assert( fixed_max_iter >= 0 );

	// object that is used to evalutate objective and constraints
	SmartPtr<CppAD::mixed::ipopt_fixed> fixed_nlp =
	new CppAD::mixed::ipopt_fixed(
		random_options,
		fixed_tolerance,
		fixed_lower,
		fixed_upper,
		fix_constraint_lower,
		fix_constraint_upper,
		fixed_in,
		random_lower,
		random_upper,
		random_in,
		mixed_object
	);

# ifndef NDEBUG
	// check derivative calculation
	bool   trace         = false;
	double relative_tol  = 1e-3;
	ok = fixed_nlp->check_grad_f(trace, relative_tol);
	if( ! ok )
	{	fatal_error("optimize_fixed: check_grad_f failed");
	}
# endif

	// Set values used for minus and plus infinity
	app->Options()->SetNumericValue(
		"nlp_lower_bound_inf", fixed_nlp->nlp_lower_bound_inf()
	);

	// variable to hold status values returned by app
	Ipopt::ApplicationReturnStatus status;

	// initialize app
	status = app->Initialize();
	ok    &= status == Ipopt::Solve_Succeeded;
	if( ! ok )
	{	fatal_error("optimize_fixed: initalization failed");
	}

	// solve the problem
	status = app->OptimizeTNLP(fixed_nlp);

	// special case where we are not trying to solve
	if( fixed_max_iter > 0 )
	{
		if( status != Ipopt::Solve_Succeeded )
		{	warning("optimize_fixed: ipopt failed to converge");
		}
		if( ! fixed_nlp->finalize_solution_ok_ )
		{	warning("optimize_fixed: solution check failed");
		}
	}
	else
	{	ok  = status == Ipopt::Maximum_Iterations_Exceeded;
		ok != status == Ipopt::Solve_Succeeded;
		if( ! ok )
		{	warning("optimize_fixed: unexpected error during zero iterations");
		}
	}
	//
	// return the entire solution including the lagrange multipliers
	CppAD::mixed::fixed_solution solution = fixed_nlp->solution();
	//
	// special case were we return fixed_in
	if( fixed_max_iter == 0 )
		solution.fixed_opt = fixed_in;
	//
	return solution;
}

