// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-25 Bradley M. Bell
// ----------------------------------------------------------------------------
# ifndef CPPAD_MIXED_CPPAD_MIXED_HPP
# define CPPAD_MIXED_CPPAD_MIXED_HPP

# include <map>
# include <cppad/cppad.hpp>
# include <cppad/mixed/sparse_hes_rcv.hpp>
# include <cppad/mixed/sparse_jac_rcv.hpp>
# include <cppad/mixed/ldlt_eigen.hpp>
# include <cppad/mixed/sparse_mat_info.hpp>
# include <cppad/mixed/configure.hpp>
# include <cppad/mixed/ldlt_cholmod.hpp>
# include <cppad/mixed/fixed_solution.hpp>
# include <cppad/mixed/typedef.hpp>
# include <cppad/mixed/warm_start_struct.hpp>
# include <cppad/mixed/trace_struct.hpp>

/*
Examples and tests that use cppad_mixed private information.
This list can be generated using the following command:
grep -P '^\t*friend *bool *::' cppad_mixed.hpp | \
   sed -e 's|^\t*friend *bool *::|extern bool |' | sort
*/
extern bool delta_ran_obj(void);
extern bool der_var_hes(void);
extern bool fix_con_eval_xam(void);
extern bool fix_con_hes_xam(void);
extern bool fix_con_jac_xam(void);
extern bool fix_like_eval_xam(void);
extern bool fix_like_hes_xam(void);
extern bool fix_like_jac_xam(void);
extern bool hes_cross_xam(void);
extern bool laplace_obj_tst(void);
extern bool laplace_obj_fun(void);
extern bool laplace_obj_hes(void);
extern bool laplace_obj_hes_xam(void);
extern bool logdet_jac_xam(void);
extern bool order2random_xam(void);
extern bool ran_con_eval_xam(void);
extern bool ran_con_jac_xam(void);
extern bool ran_hes_fun_xam(void);
extern bool ran_jac_fun_xam(void);
extern bool ran_like_hes_xam(void);
extern bool ran_obj_eval_xam(void);
extern bool ran_obj_jac_xam(void);
extern bool ran_obj_tst(void);
extern bool sample_fixed_1(void);
extern bool update_factor_xam(void);


namespace CppAD { namespace mixed {
   class optimize_random_ipopt;
   class ipopt_fixed;
   class ipopt_random;
} }

class cppad_mixed {
   friend class CppAD::mixed::optimize_random_ipopt;
   friend class CppAD::mixed::ipopt_fixed;
   friend class CppAD::mixed::ipopt_random;
   friend bool ::ran_obj_tst(void);
   friend bool ::order2random_xam(void);
   friend bool ::ran_like_hes_xam(void);
// public member functions and variables
# include <cppad/mixed/base_class.hpp>

// private member functions and variables
# include <cppad/mixed/private_base_class.hpp>
};


# include <cppad/mixed/pack.hpp>
# include <cppad/mixed/unpack.hpp>

# endif
