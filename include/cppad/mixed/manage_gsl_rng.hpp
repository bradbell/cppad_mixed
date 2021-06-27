/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_MANAGE_GSL_RNG_HPP
# define CPPAD_MIXED_MANAGE_GSL_RNG_HPP

# include <gsl/gsl_rng.h>

namespace CppAD { namespace mixed {
	size_t   new_gsl_rng(size_t seed);
	gsl_rng* get_gsl_rng(void);
	void     free_gsl_rng(void);
} }

# endif
