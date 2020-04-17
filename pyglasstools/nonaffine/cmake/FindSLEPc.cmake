## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#############################################################
# Try to find SLEPc                                         #
#                                                           #
# Once done this will define:                               #
#  SLEPC_FOUND     - system has SLEPc                       #
#  SLEPC_DIR       - SLEPc directory                        #
#  SLEPC_INC       - SLEPc include directory                #
#  SLEPC_LIB       - SLEPc library (static or dynamic)      #
#                                                           #
# Usage:                                                    #
#  find_package(SLEPc)                                      #
#                                                           #
# Setting these changes the behavior of the search          #
#  SLEPC_DIR       - SLEPc directory                        #
#############################################################

## Try to set SLEPC_DIR ##
##########################
if(NOT DEFINED SLEPC_DIR)
  set(SLEPC_DIR $ENV{SLEPC_DIR})
endif()

## Includes ##
##############
if(EXISTS "${SLEPC_DIR}/include" AND
   EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/include")
 set(SLEPC_INC "${SLEPC_DIR}/include" "${SLEPC_DIR}/${PETSC_ARCH}/include")
else()
  message(SEND_ERROR "SLEPc includes not found")
endif()

## Library ##
#############
if(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
  set(SLEPC_LIB "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.so")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
  set(SLEPC_LIB "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.a")
elseif(EXISTS "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib")
  set(SLEPC_LIB "${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib")
else()
  message(SEND_ERROR "SLEPc library not found")
endif()

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc
  "SLEPc could not be found: be sure to set SLEPC_DIR in your environment variables"
  SLEPC_LIB SLEPC_INC SLEPC_DIR)
