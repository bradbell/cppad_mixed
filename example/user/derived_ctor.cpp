// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-20 University of Washington
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
$srcthisfile%0%// BEGIN C++%// END C++%1%$$
$$

$end
*/
// BEGIN C++
# include <cppad/cppad.hpp>
# include <cppad/mixed/cppad_mixed.hpp>


namespace {
	using CppAD::log;
	using CppAD::AD;
	//
	using CppAD::mixed::d_sparse_rcv;
	using CppAD::mixed::d_vector;
	//
	class mixed_derived : public cppad_mixed {
	public:
		CppAD::vector<std::string> warning_message_;
		//
		// constructor
		mixed_derived(
			size_t                 n_fixed        ,
			size_t                 n_random       ,
			bool                   quasi_fixed    ,
			bool                   bool_sparsity  ,
			const d_sparse_rcv&    A_rcv          )
			:
			cppad_mixed(
				n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
			)
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
	size_t n_fixed     = 1;
	size_t n_random    = 0;
	bool quasi_fixed   = true;
	bool bool_sparsity = true;
	d_sparse_rcv A_rcv; // empty matrix
	//
	d_vector fixed_vec(n_fixed), random_vec(n_random);
	fixed_vec[0] = 0.0;
	//
	mixed_derived mixed_object(
		n_fixed, n_random, quasi_fixed, bool_sparsity, A_rcv
	);
	// One normaly does not store the return value of initialize
	// (size_map is included here to show how it can be displayed).
	std::map<std::string, size_t> size_map;
	size_map = mixed_object.initialize(fixed_vec, random_vec);
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
	// The following code can be used to display size_map
	// std::map<std::string, size_t>::iterator itr;
	// for(itr = size_map.begin(); itr != size_map.end(); itr++)
	//	std::cout << itr->first << " = " << itr->second << "\n";
	//
	return ok;
}
// END C++
