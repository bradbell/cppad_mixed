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
# 2DO: move cholmod_xam to test
IF( suitesparse_prefix )
	# no pkg-config for SuiteSparse so explicitly set corresponding values
	SET(suitesparse_LIBRARY_DIRS ${suitesparse_prefix}/lib )
	SET(suitesparse_LIBRARIES
		"cholmod;amd;camd;colamd;ccolamd;suitesparseconfig"
	)
	SET(suitesparse_EXAMPLES "cholmod_xam.cpp")
ELSE( suitesparse_prefix )
	SET(suitesparse_LIBRARY_DIRS "")
	SET(suitesparse_LIBRARIES    "")
	SET(suitesparse_EXAMPLES     "")
ENDIF( suitesparse_prefix )

LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
	${suitesparse_LIBRARY_DIRS}
)
# -----------------------------------------------------------------------------
# example_eigen
ADD_LIBRARY(example_eigen
	private/eigen_xam.cpp
)
SET_TARGET_PROPERTIES(example_eigen
	PROPERTIES COMPILE_FLAGS "${eigen_cxx_flags}"
)
# -----------------------------------------------------------------------------
# 2DO: move ipopt_xam_run to test.
ADD_EXECUTABLE(example EXCLUDE_FROM_ALL
	example.cpp
	${suitesparse_EXAMPLES}
	ipopt_xam_run.cpp
	private/fix_con_eval_xam.cpp
	private/fix_con_hes_xam.cpp
	private/fix_con_jac_xam.cpp
	private/fix_likelihood_eval_xam.cpp
	private/fix_likelihood_hes_xam.cpp
	private/fix_likelihood_jac_xam.cpp
	private/hes_cross_xam.cpp
	private/hes_ran_fun_xam.cpp
	private/logdet_grad_xam.cpp
	private/newton_step_xam.cpp
	private/ran_likelihood_grad_xam.cpp
	private/ran_obj_eval_xam.cpp
	private/ran_obj_grad_xam.cpp
	private/ran_obj_hes_xam.cpp
	user/data_mismatch_xam.cpp
	user/derived_xam.cpp
	user/fix_constraint_xam.cpp
	user/manage_gsl_rng_xam.cpp
	user/no_random_xam.cpp
	user/optimize_fixed_xam.cpp
	user/optimize_random_xam.cpp
)

LINK_DIRECTORIES(
	${ipopt_LIBRARY_DIRS}
)
SET_TARGET_PROPERTIES(
	example PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
TARGET_LINK_LIBRARIES(example
	example_eigen
	cppad_mixed
	cppad_mixed_eigen
	cppad_mixed
	${ipopt_LIBRARIES}
	${suitesparse_LIBRARIES}
	gsl
	gslcblas
)
ADD_CUSTOM_TARGET(check_example example DEPENDS example )
