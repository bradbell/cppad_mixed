// SPDX-License-Identifier: AGPL-3.0-or-later
// SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
// SPDX-FileContributor: 2014-22 Bradley M. Bell
// ----------------------------------------------------------------------------
# include <cassert>
# include <cppad/mixed/ipopt_app_status.hpp>

namespace CppAD { namespace mixed { // BEGIN_CPPAD_MIXED_NAMESPACE
/*
{xrst_begin ipopt_app_status dev}
{xrst_spell
   str
}

Map Ipopt Application Return Status to a String
###############################################

Syntax
******

   *str* = ``CppAD::mixed::ipopt_app_status`` ( *status* )

Prototype
*********
{xrst_spell_off}
{xrst_code cpp} */
std::string ipopt_app_status( Ipopt::ApplicationReturnStatus status )
/* {xrst_code}
{xrst_spell_on}

{xrst_end ipopt_app_status}
-----------------------------------------------------------------------------
*/
{  std::string str;
   switch(status)
   {
      case Ipopt::Solve_Succeeded:
      str = "Solve_Succeeded";
      break;

      case Ipopt::Solved_To_Acceptable_Level:
      str = "Solved_To_Acceptable_Level";
      break;

      case Ipopt::Infeasible_Problem_Detected:
      str = "Infeasible_Problem_Detected";
      break;

      case Ipopt::Search_Direction_Becomes_Too_Small:
      str = "Search_Direction_Becomes_Too_Small";
      break;

      case Ipopt::Diverging_Iterates:
      str = "Diverging_Iterates";
      break;

      case Ipopt::User_Requested_Stop:
      str = "User_Requested_Stop";
      break;

      case Ipopt::Feasible_Point_Found:
      str = "Feasible_Point_Found";
      break;

      case Ipopt::Maximum_Iterations_Exceeded:
      str = "Maximum_Iterations_Exceeded";
      break;

      case Ipopt::Restoration_Failed:
      str = "Restoration_Failed";
      break;

      case Ipopt::Error_In_Step_Computation:
      str = "Error_In_Step_Computation";
      break;

      case Ipopt::Maximum_CpuTime_Exceeded:
      str = "Maximum_CpuTime_Exceeded";
      break;

      case Ipopt::Not_Enough_Degrees_Of_Freedom:
      str = "Not_Enough_Degrees_Of_Freedom";
      break;

      case Ipopt::Invalid_Problem_Definition:
      str = "Invalid_Problem_Definition";
      break;

      case Ipopt::Invalid_Option:
      str = "Invalid_Option";
      break;

      case Ipopt::Invalid_Number_Detected:
      str = "Invalid_Number_Detected";
      break;

      case Ipopt::Unrecoverable_Exception:
      str = "Unrecoverable_Exception";
      break;

      case Ipopt::NonIpopt_Exception_Thrown:
      str = "NonIpopt_Exception_Thrown";
      break;

      case Ipopt::Insufficient_Memory:
      str = "Insufficient_Memory";
      break;

      case Ipopt::Internal_Error:
      str = "Internal_Error";
      break;

      default:
      assert(false);
      break;

   }
   return str;
}

} } // END_CPPAD_MIXED_NAMESPACE
