// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_ORDER2RANDOM_HPP
# define CPPAD_MIXED_ORDER2RANDOM_HPP

# include <cppad/mixed/cppad_mixed.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

a1_vector order2random(
   size_t                               n_fixed         ,
   size_t                               n_random        ,
   CppAD::ADFun<a1_double, double>&     jac_a1fun       ,
   const ldlt_eigen<a1_double>&         a1_ldlt_ran_hes ,
   const a1_vector&                     beta            ,
   const a1_vector&                     theta_u
);


} } // END_CPPAD_MIXED_NAMESPACE

# endif
