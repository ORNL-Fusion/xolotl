# This CMake file was taken, then modified, from:
#
# Try to find EAVL headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(EAVL)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  EAVL_PREFIX         Set this variable to the root installation of
#                      libEAVL if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  EAVL_FOUND              System has EAVL libraries and headers
#  EAVL_LIBRARIES          The EAVL library
#  EAVL_INCLUDE_DIRS       The location of EAVL headers
#  EAVL_CONFIG_DIRS       The location of EAVL config headers


find_path(EAVL_PREFIX src/common/eavl.h
    NAMES HINTS ENV EAVL_PREFIX
)

find_library(EAVL_LIBRARIES
    NAMES libeavl.a
    HINTS ${EAVL_PREFIX}/lib
)

find_path(EAVL_INCLUDE_DIRS
    NAMES common/eavl.h
    HINTS ${EAVL_PREFIX}/src
)

find_path(EAVL_CONFIG_DIRS
    NAMES eavlConfig.h
    HINTS ${EAVL_PREFIX}/config
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EAVL DEFAULT_MSG
    EAVL_LIBRARIES
    EAVL_INCLUDE_DIRS 
    EAVL_CONFIG_DIRS 
)

mark_as_advanced(
    EAVL_PREFIX_DIRS
    EAVL_LIBRARIES
    EAVL_INCLUDE_DIRS 
    EAVL_CONFIG_DIRS 
)
