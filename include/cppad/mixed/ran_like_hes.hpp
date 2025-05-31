// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_RAN_LIKE_HES_HPP
# define CPPAD_MIXED_RAN_LIKE_HES_HPP

# include <cppad/cppad.hpp>
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

a1_sparse_rcv ran_like_hes(
   size_t                               n_fixed         ,
   size_t                               n_random        ,
   CppAD::ADFun<a1_double, double>&     jac_a1fun       ,
   const sparse_rc&                     ran_hes_uu_rc   ,
   const a1_vector&                     theta_u
);


} } // END_CPPAD_MIXED_NAMESPACE

# endif
