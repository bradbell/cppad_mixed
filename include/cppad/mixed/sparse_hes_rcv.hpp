// $Id$
/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-16 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# ifndef CPPAD_MIXED_SPARSE_HES_RCV_HPP
# define CPPAD_MIXED_SPARSE_HES_RCV_HPP
/*
$begin sparse_hes_rcv$$
$spell
	rc
	rcv
	nnz
	CppAD
	const
	std
	bool
	work work
	hes
	CppAD
%$$

$section Sparse Hessian Computation Structure$$

$head Syntax$$
$codei%CppAD::mixed::sparse_hes_rcv %hes_rcv%$$

$head Private$$
This structure is an implementation detail and not part of the
$cref/CppAD::mixed/namespace/Private/$$ user API.

$head Purpose$$
This structure holds information about the Hessian
$latex \[
	H(x) = \sum_{i=0}^{m-1} w_i f_i^{(2)} (x)
\] $$
where $latex m$$ is the number of components in the function $latex f(x)$$.

$head subset$$
The field $icode%hes_rcv%.subset%$$ has prototype
$codei%
	CppAD::mixed::sparse_rcv %hes_rcv%.subset
%$$
It is empty zero when it is constructed; i.e.,
all of its sizes are zero.
After initialization it corresponds to the subset of the Hessian
that is computed.

$subhead nnz$$
We use the notation
$codei%
	%nnz% = %hes_rcv%.subset.nnz()
%$$

$subhead row$$
We use the notation
$codei%
	%row% = %hes_rcv%.subset.row()
%$$

$subhead col$$
We use the notation
$codei%
	%col% = %hes_rcv%.subset.col()
%$$

$subhead val$$
We use the notation
$codei%
	%val% = %hes_rcv%.subset.val()
%$$

$head work$$
The field $icode%hes_rcv%.work%$$ has prototype
$codei%
	CppAD::sparse_hes_work %hes_rcv%.work
%$$
It has no information when it is constructed; i.e., it is empty.
After initialization it should contain the CppAD cached information
for the call to $code sparse_hes$$ specified below.

$head Computing Sparse Hessians$$
$codei%
	%f%.sparse_hes(
		%x%,
		%w%
		%hes_rcv%.subset,
		%not_used_pattern%,
		%not_used_coloring%,
		%hes_rcv%.work
	)
%$$
Upon return,
for $icode%k% = 0 , %...%, %nnz%-1%$$,
$icode%val%[%k%]%$$ is the value of $latex H_{i,j} (x)$$
where $icode%i% = %row%[%k%]%$$ and $icode%j% = %col[%k%]%$$.

$subhead f$$
The function $icode f$$ has prototype
$codei%
	CppAD::ADFun<double> %f%
%$$

$subhead x$$
The argument $icode x$$ has prototype
$codei%
	const CppAD::vector<double>& %x%
%$$
and its size is $icode%f%.Domain()%$$.
It is the location where the Hessian is being evaluated.

$subhead w$$
The argument $icode w$$ has prototype
$codei%
	const CppAD::vector<double>& %w%
%$$
and its size is $icode%f%.Range()%$$.
It is the weighting for the components of the Hessian that is
being computed.

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
	struct sparse_hes_rcv {
		sparse_rcv                 subset;
		CppAD::sparse_hessian_work work;
	};
} }

# endif
