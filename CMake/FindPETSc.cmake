# - Try to find PETSc
# Once done this will define
#
#  PETSc_FOUND        - system has PETSc
#  PETSc_INCLUDES     - the PETSc include directories
#  PETSc_LIBRARIES    - Link these to use PETSc
#
# Setting these changes the behavior of the search
#  PETSC_DIR - directory in which PETSc resides
#  PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

find_library(PETSC_LIBRARY
    NAMES petsc
    PATHS
    ${PETSc_DIR}/lib
    ${PETSC_DIR}/lib
    $ENV{PETSC_DIR}/lib
    ${PETSC_DIR}/${PETSC_ARCH}/lib
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib
    NO_DEFAULT_PATH
)
find_library(PETSC_LIBRARY NAMES petsc)

find_path(PETSC_INCLUDE_DIR
    NAMES petsc.h
    PATHS
    ${PETSc_DIR}/include
    ${PETSC_DIR}/include
    $ENV{PETSC_DIR}/include
)

find_path(PETSCCONF_INCLUDE_DIR
    NAMES petscconf.h
    PATHS
    ${PETSc_DIR}/include
    ${PETSC_DIR}/include
    $ENV{PETSC_DIR}/include
    ${PETSC_DIR}/${PETSC_ARCH}/include
    $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include
)

get_filename_component(PETSC_LIB_DIR ${PETSC_LIBRARY} DIRECTORY)
get_filename_component(PETSC_DIR ${PETSC_LIB_DIR} DIRECTORY)
set(PETSc_DIR "${PETSC_DIR}" CACHE PATH "" FORCE)
set(PETSc_LIBRARIES "${PETSC_LIBRARY}" CACHE FILEPATH "" FORCE)
list(APPEND PETSC_INCLUDES "${PETSC_INCLUDE_DIR}" "${PETSCCONF_INCLUDE_DIR}")
list(REMOVE_DUPLICATES PETSC_INCLUDES)
set(PETSc_INCLUDES ${PETSC_INCLUDES} CACHE PATH "" FORCE)
unset(PETSC_LIBRARY CACHE)
unset(PETSC_INCLUDE_DIR CACHE)
unset(PETSCCONF_INCLUDE_DIR CACHE)

mark_as_advanced(PETSc_INCLUDES PETSc_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSc_LIBRARIES PETSc_INCLUDES
    FAIL_MESSAGE
    "PETSc could not be found. Be sure to set PETSC_DIR and PETSC_ARCH."
)

# Create interface library
add_library(PETSc UNKNOWN IMPORTED)
set_target_properties(PETSc PROPERTIES
    IMPORTED_LOCATION "${PETSc_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${PETSc_INCLUDES}"
)
