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
$codei%CppAD::mixed exception(%where%, %what%) %e%
%$$
$icode%place% = %e%.where()
%$$
$icode%description% = %e%.what()
%$$

$head Private$$
This class is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head where$$
This argument has prototype
$codei%
	const std::string& %where%
%$$
and is the name of the $code cppad_mixed$$
routine in which the exception occurred.

$head what$$
This argument has prototype
$codei%
	const std::string& %what%
%$$
and is a brief description of the exception.

$head place$$
This return value has prototype
$codei%
	std::string %place%
%$$
and is the value of $icode where$$ when $icode e$$ was constructed.

$head description$$
This argument has prototype
$codei%
	std::string %description%
%$$
and is the value of $icode what$$ when $icode e$$ was constructed.

$end
*/
# include <string>

namespace CppAD { namespace mixed {
	class exception {
	private:
		std::string where_;
		std::string what_;
	public:
		exception(const std::string& where, const std::string& what)
		: where_( where ) , what_(what)
		{ }
		std::string where(void) const
		{	return where_; }
		const std::string what(void) const
		{	return what_; }
	};
} }

# endif
