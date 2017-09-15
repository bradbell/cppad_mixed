// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_EXCEPTION_HPP
# define CPPAD_MIXED_EXCEPTION_HPP

/*
$begin exception$$
$spell
	cppad
	CppAD
	const
	std
$$

$section CppAD Mixed Exceptions$$

$head Syntax$$
$codei%CPPAD_MIXED_CALL_FATAL_ERROR
%$$
$codei%CppAD::mixed exception(%thrower%, %brief%) %e%
%$$
$icode%description% = %e%.message(%catcher%)
%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head CPPAD_MIXED_CALL_FATAL_ERROR$$
This preprocessor symbol is defined when you include this file;
see $cref/call_fatal_error/run_cmake.sh/call_fatal_error/$$

$head thrower$$
This argument has prototype
$codei%
	const std::string& %thrower%
%$$
and is the name of the routine in which the exception occurred
(the routine that threw the exception).

$head brief$$
This argument has prototype
$codei%
	const std::string& %brief%
%$$
and is a brief description of the exception.

$head catcher$$
This argument has prototype
$codei%
	const std::string& %catcher%
%$$
and is the name of the routine that caught the exception.

$head description$$
This return has prototype
$codei%
	std::string %description%
%$$
it is a message that includes
$icode catcher$$, $icode thrower$$ and $icode brief$$.

$end
*/
# include <string>

# include <cppad/mixed/configure.hpp>
# ifndef CPPAD_MIXED_CALL_FATAL_ERROR
call_fatal_error_did_not_get_defined_by_including_configure_hpp
# endif

namespace CppAD { namespace mixed {
	class exception {
	private:
		std::string thrower_;
		std::string brief_;
	public:
		exception(const std::string& thrower, const std::string& brief)
		: thrower_( thrower ) , brief_(brief)
		{ }
		const std::string message(const std::string& catcher) const
		{	std::string msg = catcher + ": " + thrower_ + ":\n" + brief_;
			return msg;
		}
	};
} }

# endif
