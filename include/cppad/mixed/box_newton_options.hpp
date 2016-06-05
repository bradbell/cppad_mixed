// $Id:$
# ifndef CPPAD_MIXED_BOX_NEWTON_OPTIONS_HPP
# define CPPAD_MIXED_BOX_NEWTON_OPTIONS_HPP
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
namespace CppAD { namespace mixed {
	// BEGIN OPTION
	struct box_newton_options {
		double tolerance;
		double direction_ratio;
		double line_ratio;
		size_t max_iter;
		size_t max_line;
		size_t print_level;
		box_newton_options(void) : // set default values
		tolerance(1e-8)       ,
		direction_ratio(0.1)  ,
		line_ratio(0.05)      ,
		max_iter(50)          ,
		max_line(10)          ,
		print_level(0)
		{}
	};
	// END OPTION

} }

# endif
