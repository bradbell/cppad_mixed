// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/ipopt/solve.hpp>
# include <cppad/mixed/cppad_mixed.hpp>
# include <cppad/mixed/box_newton.hpp>
/*
$begin optimize_random$$
$spell
	CppAD
	cppad
	vec
	const
	CppAD
	std
	Ipopt
$$

$section Optimize Random Effects$$

$head Syntax$$
$icode%random_out% =%$$
$icode%mixed_object%.optimize_random(
	%options%, %fixed_vec%, %random_lower%, %random_upper%, %random_in%
)%$$

$head Public$$
This $code cppad_mixed$$ member function is $cref public$$.

$head Purpose$$
This routine maximizes the
$cref/random likelihood/ran_likelihood/$$
corresponding to the object $icode mixed_object$$.

$head mixed_object$$
We use $cref/mixed_object/derived_ctor/mixed_object/$$
to denote an object of a class that is
derived from the $code cppad_mixed$$ base class.

$head options$$
The argument $icode options$$ has prototype
$codei%
	const std::string& %options%
%$$
and is the $cref ipopt_options$$ for optimizing the random effects.

$head fixed_vec$$
This argument has prototype
$codei%
	const CppAD::vector<double>& %fixed_vec%
%$$
It specifies the value of the
$cref/fixed effects/cppad_mixed/Notation/Fixed Effects, theta/$$
vector $latex \theta$$.

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

$head random_out$$
The return value  has prototype
$codei%
	CppAD::vector<double> %random_out%
%$$
It is the final value (obtained by optimization) of the
$cref/random effects/cppad_mixed/Notation/Random Effects, u/$$
vector $latex u$$.

$children%
	example/user/optimize_random_xam.cpp
%$$
$head Example$$
The file $cref optimize_random_xam.cpp$$ contains an example
and test of this procedure.
It returns true, if the test passes, and false otherwise.

$end
*/


namespace CppAD{ namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

// ----------------------------------------------------------------------------
// helper class used by optimize_random
class optimize_random_box_newton {
	typedef CppAD::vector<double> d_vector;
private:
	const size_t         n_random_;     // number of random effects
	const size_t         n_fixed_;      // number of fixed effects
	const d_vector       fixed_vec_;    // value of fixed effects
	size_t               n_iter_;       // number of iterations so far
	d_vector             random_vec_;   // previous random_vec
	cppad_mixed&         mixed_object_; // cppad_mixed for this problem
public:
	// constructor
	optimize_random_box_newton(
		size_t            n_random      ,
		const d_vector&   fixed_vec     ,
		cppad_mixed&      mixed_object  ) :
	n_random_( n_random )          ,
	n_fixed_( fixed_vec.size() )   ,
	fixed_vec_( fixed_vec )        ,
	random_vec_( n_random )         ,
	mixed_object_( mixed_object )
	{	n_iter_ = 0; }
	// fun
	double fun(const d_vector& random_vec)
	{	assert( random_vec.size() == n_random_ );
		d_vector both(n_fixed_ + n_random_);
		mixed_object_.pack(fixed_vec_, random_vec, both);
		// f(theta, u)
		d_vector vec = mixed_object_.ran_like_fun_.Forward(0, both);
		assert( vec.size() == 1 );
		random_vec_ = random_vec;
		return vec[0];
	}
	// grad
	d_vector grad(const d_vector& random_vec)
	{	assert( random_vec.size() == n_random_ );
# ifndef NDEBUG
		for(size_t i = 0; i < n_random_; i++)
			assert( random_vec[i] == random_vec_[i] );
# endif
		// use fact that previous call to fun was with same random_vec
		d_vector w(1);
		w[0] = 1.0;
		d_vector grad_both = mixed_object_.ran_like_fun_.Reverse(1, w);
		assert( grad_both.size() == n_fixed_ + n_random_ );
		d_vector grad_fixed(n_fixed_), grad_random(n_random_);
		mixed_object_.unpack(grad_fixed, grad_random, grad_both);
		return grad_random;
	}
	// solve
	d_vector solve(
		const d_vector& random_vec  ,
		const d_vector& p           )
	{
# ifndef NDEBUG
		for(size_t i = 0; i < n_random_; i++)
			assert( random_vec[i] == random_vec_[i] );
# endif
		// update the Hessian factor mixed_object_.ldlt_ran_hes_
		// once every 10 iterations
		if( n_iter_ % 5 == 0 )
			mixed_object_.update_factor(fixed_vec_, random_vec);
		n_iter_++;
		// solve H * d = p
		d_vector d(n_random_);
		CppAD::vector<size_t> row(n_random_);
		for(size_t i = 0; i < n_random_; i++)
			row[i] = i;
		mixed_object_.ldlt_ran_hes_.solve_H(row, p, d);
		return d;
	}
};

} } // END_CPPAD_MIXED_NAMESPACE

// ----------------------------------------------------------------------------
// optimize_random
CppAD::vector<double> cppad_mixed::optimize_random(
	const std::string& not_used , // will be changed to box_newton_options
	const d_vector&    fixed_vec       ,
	const d_vector&    random_lower    ,
	const d_vector&    random_upper    ,
	const d_vector&    random_in       )
{
	// make sure initialize has been called
	if( ! initialize_done_ )
	{	std::string error_message =
		"cppad_mixed::initialize was not called before optimize_random";
		fatal_error(error_message);
	}
	if( ! init_ran_like_done_ )
	{	std::string error_message =
		"cppad_mixed::optimize_random there are no random effects";
		fatal_error(error_message);
	}

	// number of fixed and random effects
	assert( n_fixed_  == fixed_vec.size() );
	assert( n_random_ == random_in.size() );
	assert( n_random_ == random_lower.size() );
	assert( n_random_ == random_upper.size() );

	// call back object used by optimizer
	CppAD::mixed::optimize_random_box_newton objective(
		n_random_, fixed_vec, *this
	);
	//
	// call optimizer
	CppAD::mixed::box_newton_options options;
	options.print_level = 0;
	options.tolerance   = 1e-10;
	d_vector random_out(n_random_);
	CppAD::mixed::box_newton_status status = CppAD::mixed::box_newton(
		options, objective, random_lower, random_upper, random_in, random_out
	);
	switch( status )
	{	case CppAD::mixed::box_newton_ok_enum:
		break;

		case CppAD::mixed::box_newton_max_iter_enum:
		fatal_error("optmize_random: box_newton_max_iter");
		break;

		case CppAD::mixed::box_newton_max_line_enum:
		fatal_error("optmize_random: box_newton_max_line");
		break;

		default:
		assert(false);
		break;
	}

	return random_out;
}



