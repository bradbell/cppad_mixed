# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
#           Copyright (C) 2014-16 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# Build C++ Examples / Tests
#
LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
	${suitesparse_LIBRARY_DIRS}
)
# -----------------------------------------------------------------------------
# 2DO: move ipopt_xam_run to test.
# 2DO: move cholmod_xam to test
ADD_EXECUTABLE(example EXCLUDE_FROM_ALL
	example.cpp
	private/cholmod_solve_xam.cpp
	private/cholmod_solve2_xam.cpp
	private/cholmod_xam.cpp
	ipopt_xam_run.cpp
	private/eigen_xam.cpp
	private/fix_con_eval_xam.cpp
	private/fix_con_hes_xam.cpp
	private/fix_con_jac_xam.cpp
	private/fix_like_eval_xam.cpp
	private/fix_like_hes_xam.cpp
	private/fix_like_jac_xam.cpp
	private/hes_cross_xam.cpp
	private/ran_hes_fun_xam.cpp
	private/logdet_jac_xam.cpp
	private/newton_step_xam.cpp
	private/ran_con_eval_xam.cpp
	private/ran_con_jac_xam.cpp
	private/ran_like_jac_xam.cpp
	private/ran_obj_eval_xam.cpp
	private/ran_obj_jac_xam.cpp
	private/ran_objcon_hes_xam.cpp
	private/undetermined_xam.cpp
	private/update_factor_xam.cpp
	user/abs_density_xam.cpp
	user/data_mismatch_xam.cpp
	user/derived_ctor_xam.cpp
	user/fix_constraint_xam.cpp
	user/fix_likelihood_xam.cpp
	user/information_mat_xam.cpp
	user/lasso_xam.cpp
	user/manage_gsl_rng_xam.cpp
	user/no_random_xam.cpp
	user/optimize_fixed_xam.cpp
	user/optimize_random_xam.cpp
	user/ran_constraint_xam.cpp
	user/ran_likelihood_xam.cpp
	user/ran_likelihood_hes_xam.cpp
	user/ran_likelihood_jac_xam.cpp
	user/sample_fixed_xam.cpp
)

LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
)
SET_TARGET_PROPERTIES(
	example PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(example
	cppad_mixed
	${ipopt_LIBRARIES}
	${suitesparse_LIBRARIES}
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(check_example example DEPENDS example )
