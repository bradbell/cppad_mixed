/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_JAC_RCV_HPP
# define CPPAD_MIXED_SPARSE_JAC_RCV_HPP
/*
$begin sparse_jac_rcv$$
$spell
	jac
	nnz
	rc
	rcv
	CppAD
	const
	std
	bool
	work work
	hes
	CppAD
	Jacobian
	Jacobians
%$$

$section Sparse Jacobian Computation Structure$$

$head Syntax$$
$codei%CppAD::mixed::sparse_jac_rcv %jac_rcv%$$

$head Private$$
This structure is an implementation detail and not part of the
CppAD Mixed user API.

$head Purpose$$
This structure holds information about a specific Jacobian
$latex \[
	J(x) = f^{(1)} (x)
\] $$

$head subset$$
The field $icode%jac_rcv%.subset%$$ has prototype
$codei%
	CppAD::mixed::d_sparse_rcv %jac_rcv%.subset
%$$
It is empty zero when it is constructed; i.e.,
all of its sizes are zero.
After initialization it corresponds to the subset of the Jacobian
that is computed.

$subhead nnz$$
We use the notation
$codei%
	%nnz% = %jac_rcv%.subset.nnz()
%$$

$subhead row$$
We use the notation
$codei%
	%row% = %jac_rcv%.subset.row()
%$$

$subhead col$$
We use the notation
$codei%
	%col% = %jac_rcv%.subset.col()
%$$

$subhead val$$
We use the notation
$codei%
	%val% = %jac_rcv%.subset.val()
%$$

$head forward$$
The field $icode%jac_rcv%.forward%$$ has prototype
$codei%
	bool %jac_rcv%.forward%
%$$

$head work$$
The field $icode%jac_rcv%.work%$$ has prototype
$codei%
	CppAD::sparse_jac_work %jac_rcv%.work
%$$
It has no information when it is constructed; i.e., it is empty.
After initialization it should contain the CppAD cached information
for the call to $code sparse_jac_for$$ or $code sparse_jac_rev$$
specified below.

$head Computing Sparse Jacobians$$
Upon return (from $code sparse_jac_for$$ or $code sparse_jac_rev$$),
for $icode%k% = 0 , %...%, %nnz%-1%$$,
$icode%val%[%k%]%$$ is the value of $latex J_{i,j} (x)$$
where $icode%i% = %row%[%k%]%$$ and $icode%j% = %col[%k%]%$$.

$subhead Forward Mode$$
If $icode%jac_rcv%.forward%$$ is true,
$codei%
	%f%.sparse_jac_for(
		%group_max%,
		%x%,
		%jac_rcv%.subset,
		%not_used_pattern%,
		%not_used_coloring%,
		%jac_rcv%.work
	)
%$$
computes the Jacobian values.

$subhead Reverse Mode$$
If $icode%jac_rcv%.forward%$$ is false,
$codei%
	%f%.sparse_jac_for(
		%x%,
		%jac_rcv%.subset,
		%not_used_pattern%,
		%not_used_coloring%,
		%jac_rcv%.work
	)
%$$
computes the Jacobian values.

$subhead f$$
The function $icode f$$ has prototype
$codei%
	CppAD::ADFun<double> %f%
%$$

$subhead group_max$$
The parameter $icode group_max$$ has prototype
$codei%
	size_t %group_max%
%$$
and must be greater than zero.
It specifies the maximum number of directions to group during
a single forward sweep.
This uses separate memory for each direction (more memory),
but my be significantly faster.

$subhead x$$
The argument $icode x$$ has prototype
$codei%
	const CppAD::vector<double>& %x%
%$$
and its size is $icode%f%.Domain()%$$.
It is the location where the Jacobian is being evaluated.

$subhead not_used_pattern$$
This argument has the following prototype
$codei%
	const CppAD::mixed::sparse_rc& %not_used_pattern%
%$$
It is not used and hence its value does not matter.

$subhead not_used_coloring$$
This argument has the following prototype
$codei%
	const std::string& %not_used_coloring%
%$$
It is not used and hence its value does not matter.

$end
*/
# include <cppad/mixed/typedef.hpp>

namespace CppAD { namespace mixed {
	struct sparse_jac_rcv {
		bool                   forward;
		d_sparse_rcv           subset;
		CppAD::sparse_jac_work work;
	};
} }

# endif
