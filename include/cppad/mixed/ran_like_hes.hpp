// $Id:$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-18 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_RAN_LIKE_HES_HPP
# define CPPAD_MIXED_RAN_LIKE_HES_HPP

# include <cppad/cppad.hpp>
# include <cppad/mixed/typedef.hpp>
# include <cppad/mixed/ldlt_eigen.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE

a1_sparse_rcv ran_like_hes(
	size_t                       n_fixed         ,
	size_t                       n_random        ,
	CppAD::ADFun<a1_double>&     jac_a1fun       ,
	const sparse_rc&             ran_hes_uu_rc   ,
	const a1_vector&             theta_u
);


} } // END_CPPAD_MIXED_NAMESPACE

# endif
