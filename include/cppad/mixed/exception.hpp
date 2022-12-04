// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_EXCEPTION_HPP
# define CPPAD_MIXED_EXCEPTION_HPP

/*
{xrst_begin exception}

CppAD Mixed Exceptions
######################

Syntax
******

| ``CppAD::mixed exception`` ( *thrower* , *brief* ) *e*
| *description* = *e* . ``message`` ( *catcher* )

thrower
*******
This argument has prototype

   ``const std::string&`` *thrower*

and is the name of the routine in which the exception occurred
(the routine that threw the exception).

brief
*****
This argument has prototype

   ``const std::string&`` *brief*

and is a brief description of the exception.

catcher
*******
This argument has prototype

   ``const std::string&`` *catcher*

and is the name of the routine that caught the exception.

description
***********
This return has prototype

   ``std::string`` *description*

it is a message that includes
*catcher* , *thrower* and *brief* .

{xrst_end exception}
*/
# include <string>
# include <cppad/mixed/configure.hpp>

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
      {  std::string msg = catcher + ": " + thrower_ + ":\n" + brief_;
         return msg;
      }
   };
} }

# endif
