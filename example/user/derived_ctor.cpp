// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
/*
$begin derived_ctor.cpp$$
$spell
	CppAD
	cppad
	interp
	xam
$$

$section mixed_cppad Derived Class: Example and Test$$

$code
$srcfile%example/user/derived_ctor.cpp
	%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>


namespace {
	using CppAD::vector;
	using CppAD::log;
	using CppAD::AD;
	using CppAD::mixed::sparse_mat_info;
	//
	typedef AD<double>    a1_double;
	typedef AD<a1_double> a2_double;

	class mixed_derived : public cppad_mixed {
	public:
		vector<std::string> warning_message_;
		//
		// constructor
		mixed_derived(
			size_t n_fixed                    ,
			size_t n_random                   ,
			bool   quasi_fixed    )
			:
			cppad_mixed(n_fixed, n_random, quasi_fixed)
		{ }
		// example changing fatal error handler
		virtual void fatal_error(const std::string& error_message)
		{	// throw a std:string
			throw error_message;
		}
		//
		// example changing warning handler
		virtual void warning(const std::string& warning_message)
		{	warning_message_.push_back( warning_message );
		}
	};
}

bool derived_ctor_xam(void)
{
	//
	size_t n_fixed   = 1;
	size_t n_random  = 0;
	bool quasi_fixed = true;
	//
	CppAD::mixed::sparse_mat_info A_info;
	vector<double> fixed_vec(n_fixed), random_vec(n_random);
	fixed_vec[0] = 0.0;
	//
	mixed_derived mixed_object(n_fixed, n_random, quasi_fixed);
	mixed_object.initialize(fixed_vec, random_vec, A_info);
	//
	bool ok = false;
	try
	{	// warnings
		mixed_object.warning("first warning");
		mixed_object.warning("second warning");

		// fatal error
		mixed_object.fatal_error("only fatal error");
	}
	catch ( std::string error_message)
	{	ok  = true;
		ok &= error_message == "only fatal error";
		//
		ok &= mixed_object.warning_message_.size() == 2;
		ok &= mixed_object.warning_message_[0] == "first warning";
		ok &= mixed_object.warning_message_[1] == "second warning";
	}
	return ok;
}
// END C++
