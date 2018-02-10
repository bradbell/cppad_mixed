# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-18 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
ADD_LIBRARY(cppad_mixed
	cholmod/constructor.cpp
	cholmod/init.cpp
	cholmod/inv.cpp
	cholmod/logdet.cpp
	cholmod/solve_H.cpp
	cholmod/sim_cov.cpp
	cholmod/update.cpp
	derived_ctor.cpp
	eigen/ldlt_eigen.cpp
	eigen/information_mat.cpp
	eigen/ldlt_inv.cpp
	eigen/newton_step.cpp
	eigen/sample_fixed.cpp
	eigen/sparse_low_tri_sol.cpp
	eigen/sparse_up_tri_sol.cpp
	sample_random.cpp
	eigen/undetermined.cpp
	fix_con_eval.cpp
	fix_con_eval.cpp
	fix_con_hes.cpp
	fix_con_jac.cpp
	fix_like_eval.cpp
	fix_like_hes.cpp
	fix_like_jac.cpp
	init_ldlt_ran_hes.cpp
	init_fix_con.cpp
	init_fix_like.cpp
	init_hes_cross.cpp
	initialize.cpp
	init_ran_jac.cpp
	init_ran_hes.cpp
	init_ran_like.cpp
	init_laplace_obj.cpp
	init_laplace_obj_hes.cpp
	ipopt_app_status.cpp
	ipopt_fixed.cpp
	ipopt_fixed/eval_grad_f.cpp
	ipopt_fixed/eval_f.cpp
	ipopt_fixed/get_starting_point.cpp
	ipopt_fixed/get_bounds_info.cpp
	ipopt_fixed/get_nlp_info.cpp
	ipopt_fixed/ctor.cpp
	ipopt_random.cpp
	logdet_jac.cpp
	manage_gsl_rng.cpp
	optimize_fixed.cpp
	optimize_random.cpp
	ran_con_eval.cpp
	ran_con_jac.cpp
	ran_like_jac.cpp
	laplace_obj_hes.cpp
	ran_obj_eval.cpp
	ran_obj_jac.cpp
	update_factor.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed
	PROPERTIES COMPILE_FLAGS "${cppad_cxx_flags}"
)
INSTALL(
	TARGETS cppad_mixed
	DESTINATION ${cppad_prefix}/${cmake_libdir}
)
