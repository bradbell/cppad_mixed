# SPDX-License-Identifier: AGPL-3.0-or-later
# SPDX-FileCopyrightText: University of Washington <https://www.washington.edu>
# SPDX-FileContributor: 2014-25 Bradley M. Bell
# ----------------------------------------------------------------------------
# add_to_set(variable_set elememt_set)
#
# ${variable_set}: (in/out)
# A list containing the values in the set. The empty set is represented
# by the value ""; Elements in the set are separated by the semi-colon.
#
# element_set: (in)
# is the element we are adding to the set. 
# If this element is "" or if it is already in the set, 
# the set will not change. Otherise this element is added at the
# end of the list that represents the set.
#
# This macro using variables that begin with add_to_set
MACRO(add_to_set variable_set element_set)
   IF( "${${variable_set}}" STREQUAL "" )
      SET( ${variable_set} "${element_set}" )
   ELSE( )
      LIST(FIND "${${variable_set}} " "${element_set}" add_to_set_index)
      IF( "${add_to_set_index}" STREQUAL "-1" )
         SET( ${variable_set} "${${variable_set}};${element_set}" )
      ENDIF( )
   ENDIF( )
ENDMACRO( )
   
 
