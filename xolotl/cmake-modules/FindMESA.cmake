# This CMake file was taken, then modified, from:
#
# Try to find MESA headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(MESA)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  MESA_PREFIX         Set this variable to the root installation of
#                      libMESA if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  MESA_FOUND              System has MESA libraries and headers
#  MESA_LIBRARIES          The MESA library
#  MESA_INCLUDE_DIRS       The location of MESA headers


find_path(MESA_PREFIX include/GL/osmesa.h
    NAMES HINTS ENV MESA_PREFIX
)

find_library(MESA_LIBRARIES
    NAMES libOSMesa.a
    HINTS ${MESA_PREFIX}/lib
)

find_path(MESA_INCLUDE_DIRS
    NAMES GL/osmesa.h
    HINTS ${MESA_PREFIX}/include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MESA DEFAULT_MSG
    MESA_LIBRARIES
    MESA_INCLUDE_DIRS 
)

mark_as_advanced(
    MESA_PREFIX_DIRS
    MESA_LIBRARIES
    MESA_INCLUDE_DIRS 
)
