// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_MANAGE_GSL_RNG_HPP
# define CPPAD_MIXED_MANAGE_GSL_RNG_HPP

# include <gsl/gsl_rng.h>

namespace CppAD { namespace mixed {
   size_t   new_gsl_rng(size_t seed);
   gsl_rng* get_gsl_rng(void);
   void     free_gsl_rng(void);
} }

# endif
