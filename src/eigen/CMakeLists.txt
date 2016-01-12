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
ADD_LIBRARY(cppad_mixed_eigen
	cholesky.cpp
	init_hes_ran.cpp
	init_ran_con.cpp
	init_ran_objcon.cpp
	logdet_jac.cpp
	newton_step.cpp
	ran_con_eval.cpp
	ran_con_jac.cpp
	ran_obj_eval.cpp
	ran_obj_jac.cpp
	update_factor.cpp
)
SET_TARGET_PROPERTIES(cppad_mixed_eigen
	PROPERTIES COMPILE_FLAGS "${extra_cxx_flags}"
)
INSTALL(
	TARGETS cppad_mixed_eigen
	DESTINATION ${cppad_prefix}/${cmake_libdir}
)
