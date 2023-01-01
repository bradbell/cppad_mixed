/* --------------------------------------------------------------------------
cppad_mixed: C++ Laplace Approximation of Mixed Effects Models
          Copyright (C) 2014-21 University of Washington
             (Bradley M. Bell bradbell@uw.edu)

This program is distributed under the terms of the
	     GNU Affero General Public License version 3.0 or later
see http://www.gnu.org/licenses/agpl.txt
-------------------------------------------------------------------------- */
# include <cppad/mixed/ipopt_fixed.hpp>
# include <coin-or/IpIpoptCalculatedQuantities.hpp>
# include <coin-or/IpOrigIpoptNLP.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/* %$$
-------------------------------------------------------------------------------
$begin ipopt_fixed_intermediate_callback$$
$spell
	Ipopt
	iter
	obj
	inf
	pr
	du
	ls
	ip_cq
	lg
	rg
	optimizer
	mu
	enum
$$

$nospell
$bold This is cppad_mixed--20220519 documentation:$$ Here is a link to its
$href%https://cppad-mixed.readthedocs.io%current documentation%$$.
$$
$section Get Optimizer Trace Information From Ipopt$$

$head Syntax$$
$icode%ok% = intermediate_callback(
	%mode%,
	%iter%,
	%obj_value%,
	%inf_pr%,
	%inf_du%,
	%mu%,
	%d_norm%,
	%regularization_size%,
	%alpha_du%,
	%alpha_pr%,
	%ls_trials%,
	%ip_data%,
	%ip_cq%
)%$$

$head mode$$
This is an $code enum$$ value and equal to
$code RegularMode$$ or $code RestorationPhaseMode$$.

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

$head ok$$
If set to false, the optimizer will terminate with
$code USER_REQUESTED_STOP$$ as the finalize_solution
$cref/status/ipopt_xam_finalize_solution/status/$$.

$head Source$$
$srccode@cpp@ */
bool ipopt_fixed::intermediate_callback(
	AlgorithmMode               mode                 ,   // in
	Index                       iter                 ,   // in
	Number                      obj_value            ,   // in
	Number                      inf_pr               ,   // in
	Number                      inf_du               ,   // in
	Number                      mu                   ,   // in
	Number                      d_norm               ,   // in
	Number                      regularization_size  ,   // in
	Number                      alpha_du             ,   // in
	Number                      alpha_pr             ,   // in
	Index                       ls_trials            ,   // in
	const IpoptData*            ip_data              ,   // in
	IpoptCalculatedQuantities*  ip_cq                )   // in
{	assert( solution_.trace_vec.size() == size_t(iter) );
	//
	bool restoration = Ipopt::GetRawPtr(ip_cq->GetIpoptNLP()) == nullptr;
	//
	trace_struct trace;
	trace.iter                = size_t(iter);
	trace.obj_value           = obj_value;
	trace.inf_pr              = inf_pr;
	trace.inf_du              = inf_du;
	trace.mu                  = mu;
	trace.d_norm              = d_norm;
	trace.regularization_size = regularization_size;
	trace.alpha_du            = alpha_du;
	trace.alpha_pr            = alpha_pr;
	trace.ls_trials           = size_t(ls_trials);
	trace.restoration         = restoration;
	//
	solution_.trace_vec.push_back(trace);
	return true;
}

} } // END_CPPAD_MIXED_NAMESPACE
