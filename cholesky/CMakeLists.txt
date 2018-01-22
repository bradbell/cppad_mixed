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
# Build C++ Examples / Tests
# -----------------------------------------------------------------------------
# bin/test_one.sh expects example.cpp on next line, last line to be just ')'
ADD_EXECUTABLE(cholesky EXCLUDE_FROM_ALL cholesky.cpp
	cholesky.cpp
	sparse_ad_cholesky.cpp
	example/sparse_ad_chol_eq.cpp
	example/sparse_ad_chol_eval.cpp
	example/sparse_ad_chol_perm.cpp
	example/sparse_ad_chol_sp1.cpp
	example/sparse_ad_chol_sp2.cpp
	example/sparse_ad_chol_var.cpp
)
SET_TARGET_PROPERTIES(
	cholesky PROPERTIES COMPILE_FLAGS "${cppad_cxx_flags}"
)
TARGET_LINK_LIBRARIES(cholesky
	cppad_mixed
	${suitesparse_library_path_list}
)
ADD_CUSTOM_TARGET(check_cholesky cholesky DEPENDS cholesky )
