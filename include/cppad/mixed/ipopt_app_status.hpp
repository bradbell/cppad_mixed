// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_IPOPT_APP_STATUS_HPP
# define CPPAD_MIXED_IPOPT_APP_STATUS_HPP

# include <string>
# include <coin-or/IpIpoptApplication.hpp>

namespace CppAD { namespace mixed {
	std::string ipopt_app_status( Ipopt::ApplicationReturnStatus status );
} }

# endif
