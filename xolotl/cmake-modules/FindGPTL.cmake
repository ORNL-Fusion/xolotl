# This CMake file was taken, then modified, from:
#
# http://git.icir.org/binpacpp.git/blob_plain/feae71bbc70f46a614f293f9db88d276a6c94379:/cmake/FindPAPI.cmake
#
# Try to find GPTL headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(GPTL)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  GPTL_PREFIX         Set this variable to the root installation of
#                      libpapi if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  GPTL_FOUND              System has GPTL libraries and headers
#  GPTL_LIBRARIES          The GPTL library
#  GPTL_INCLUDE_DIRS       The location of GPTL headers

find_path(GPTL_PREFIX include/gptl.h
    NAMES HINTS ENV GPTL_PREFIX
)

#find_path(GPTL_PREFIX include/gptl.h
#    HINTS ENV GPTL_PREFIX
#)

find_library(GPTL_LIBRARIES
    NAMES libgptl_pmpi.a gptl
#    HINTS ${GPTL_PREFIX} ${GPTL_PREFIX}/lib ${HILTIDEPS}/lib
    HINTS ${GPTL_PREFIX}/lib ${HILTIDEPS}/lib
)

find_path(GPTL_INCLUDE_DIRS
    NAMES gptl.h
    HINTS ${GPTL_PREFIX} ${GPTL_PREFIX}/include ${HILTIDEPS}/include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPTL DEFAULT_MSG
    GPTL_LIBRARIES
    GPTL_INCLUDE_DIRS
)

mark_as_advanced(
    GPTL_PREFIX_DIRS
    GPTL_LIBRARIES
    GPTL_INCLUDE_DIRS
)
