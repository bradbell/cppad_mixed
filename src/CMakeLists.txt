# $Id$
#  --------------------------------------------------------------------------
# cppad_mixed: Estimating Disease Rates as Functions of Age and Time
#           Copyright (C) 2014-15 University of Washington
#              (Bradley M. Bell bradbell@uw.edu)
#
# This program is distributed under the terms of the
#	     GNU Affero General Public License version 3.0 or later
# see http://www.gnu.org/licenses/agpl.txt
# ---------------------------------------------------------------------------
# cppad_mixed
ADD_LIBRARY(cppad_mixed
	constraint_eval.cpp
	constraint_eval.cpp
	constraint_jac.cpp
	constraint_hes.cpp
	ran_like_grad.cpp
	ranobj_hes.cpp
	initialize.cpp
	ipopt_fixed.cpp
	manage_gsl_rng.cpp
	optimize_fixed.cpp
	optimize_random.cpp
	fix_like_eval.cpp
	fix_like_jac.cpp
	fix_like_hes.cpp
	init_fix_like.cpp
	init_constraint.cpp
	init_hes_cross.cpp
	init_hes_ranobj.cpp
	init_ran_like.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed
	PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
# ---------------------------------------------------------------------------
# cppad_mixed_eigen
ADD_LIBRARY(cppad_mixed_eigen
	chol_hes_ran.cpp
	init_hes_ran.cpp
	init_ranobj.cpp
	newton_step.cpp
	ranobj_eval.cpp
	logdet_grad.cpp
	ranobj_grad.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed_eigen
	PROPERTIES COMPILE_FLAGS "${eigen_cxx_flags}"
)
# ---------------------------------------------------------------------------
INSTALL(
	TARGETS cppad_mixed cppad_mixed_eigen
		DESTINATION ${cppad_prefix}/${cppad_libdir}
)
INSTALL(
	FILES
	${cppad_mixed_SOURCE_DIR}/include/cppad/mixed/cppad_mixed.hpp
	${cppad_mixed_SOURCE_DIR}/include/cppad/mixed/manage_gsl_rng.hpp
	DESTINATION ${cppad_mixed_prefix}/include/cppad
)
