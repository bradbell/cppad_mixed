// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_IPOPT_APP_STATUS_HPP
# define CPPAD_MIXED_IPOPT_APP_STATUS_HPP

# include <string>
# include <coin/IpIpoptApplication.hpp>

namespace CppAD { namespace mixed {
	std::string ipopt_app_status( Ipopt::ApplicationReturnStatus status );
} }

# endif
