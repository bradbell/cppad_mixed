# ifndef CPPAD_MIXED_TRACE_STRUCT_HPP
# define CPPAD_MIXED_TRACE_STRUCT_HPP
// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------

/*
$begin trace_struct$$
$spell
	Ipopt
	CppAD
	struct
	iter
	obj
	inf
	pr
	du
	mu
	lg
	rg
	ls
$$

$section Ipopt Trace Information$$

$head Syntax$$
$codei%CppAD::mixed::trace_struct %trace%$$

$head Prototype$$
$srcthisfile%0%// BEGIN_PROTOTYPE%// END_PROTOTYPE%1%$$

$head iter$$
See ipopt trace $cref/iter/ipopt_trace/iter/$$.

$head obj_value$$
See ipopt trace $cref/objective/ipopt_trace/objective/$$.

$head inf_pr$$
See ipopt trace $cref/inf_pr/ipopt_trace/inf_pr/$$.

$head inf_du$$
See ipopt trace $cref/inf_du/ipopt_trace/inf_du/$$.

$head mu$$
See ipopt trace $cref/lg(mu)/ipopt_trace/lg(mu)/$$.

$head d_norm$$
See ipopt trace $cref/||d||/ipopt_trace/||d||/$$.

$head regularization_size$$
See ipopt trace $cref/lg(rg)/ipopt_trace/lg(rg)/$$.

$head alpha_du$$
See ipopt trace $cref/alpha_du/ipopt_trace/alpha_du/$$.

$head alpha_pr$$
See ipopt trace $cref/alpha_pr/ipopt_trace/alpha_pr/$$.

$head ls_trials$$
See ipopt trace $cref/ls/ipopt_trace/ls/$$.

$head restoration$$
Is ipopt currently in restoration mode.

$end
------------------------------------------------------------------------------
*/
namespace CppAD { namespace mixed {
	// BEGIN_PROTOTYPE
	struct trace_struct {
		size_t iter;
		double obj_value;
		double inf_pr;
		double inf_du;
		double mu;
		double d_norm;
		double regularization_size;
		double alpha_du;
		double alpha_pr;
		size_t ls_trials;
		bool   restoration;
	};
	// END_PROTOTYPE
} }


# endif
