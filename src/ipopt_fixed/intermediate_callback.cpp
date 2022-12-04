// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cppad/mixed/ipopt_fixed.hpp>
# include <coin-or/IpIpoptCalculatedQuantities.hpp>
# include <coin-or/IpOrigIpoptNLP.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/* %$$
-------------------------------------------------------------------------------
{xrst_begin ipopt_fixed_intermediate_callback}
{xrst_spell
   cq
   enum
   iter
   ls
   optimizer
   pr
}

Get Optimizer Trace Information From Ipopt
##########################################

Syntax
******

| *ok* = ``intermediate_callback`` (
| |tab| *mode* ,
| |tab| *iter* ,
| |tab| *obj_value* ,
| |tab| *inf_pr* ,
| |tab| *inf_du* ,
| |tab| *mu* ,
| |tab| *d_norm* ,
| |tab| *regularization_size* ,
| |tab| *alpha_du* ,
| |tab| *alpha_pr* ,
| |tab| *ls_trials* ,
| |tab| *ip_data* ,
| |tab| *ip_cq*
| )

mode
****
This is an ``enum`` value and equal to
``RegularMode`` or ``RestorationPhaseMode`` .

iter
****
See ipopt trace :ref:`ipopt_trace@iter` .

obj_value
*********
See ipopt trace :ref:`ipopt_trace@objective` .

inf_pr
******
See ipopt trace :ref:`ipopt_trace@inf_pr` .

inf_du
******
See ipopt trace :ref:`ipopt_trace@inf_du` .

mu
**
See ipopt trace :ref:`ipopt_trace@lg(mu)` .

d_norm
******
See ipopt trace :ref:`ipopt_trace@||d||` .

regularization_size
*******************
See ipopt trace :ref:`ipopt_trace@lg(rg)` .

alpha_du
********
See ipopt trace :ref:`ipopt_trace@alpha_du` .

alpha_pr
********
See ipopt trace :ref:`ipopt_trace@alpha_pr` .

ls_trials
*********
See ipopt trace :ref:`ipopt_trace@ls` .

ok
**
If set to false, the optimizer will terminate with
``USER_REQUESTED_STOP`` as the finalize_solution
:ref:`ipopt_xam_finalize_solution@status` .

Source
******
{xrst_spell_off}
{xrst_code cpp} */
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
{  assert( solution_.trace_vec.size() == size_t(iter) );
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
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_fixed_intermediate_callback}
*/

} } // END_CPPAD_MIXED_NAMESPACE
