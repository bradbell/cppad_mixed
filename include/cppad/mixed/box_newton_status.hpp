// $Id:$
# ifndef CPPAD_MIXED_BOX_NEWTON_STATUS_HPP
# define CPPAD_MIXED_BOX_NEWTON_STATUS_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
namespace CppAD { namespace mixed {
	// BEGIN STATUS
	enum box_newton_status {
		box_newton_ok_enum       , // x_out is ok
		box_newton_max_iter_enum , // max number of iterations reached
		box_newton_max_line_enum   // max number of line search steps reached
	};
	// END STATUS
} }

# endif
