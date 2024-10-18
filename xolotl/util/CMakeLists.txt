set(XOLOTL_UTIL_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_UTIL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_UTIL_HEADER_DIR ${XOLOTL_UTIL_INCLUDE_DIR}/xolotl/util)

set(XOLOTL_UTIL_HEADERS
    ${XOLOTL_UTIL_HEADER_DIR}/Array.h
    ${XOLOTL_UTIL_HEADER_DIR}/DoInOrder.h
    ${XOLOTL_UTIL_HEADER_DIR}/Filesystem.h
    ${XOLOTL_UTIL_HEADER_DIR}/IStepSequence.h
    ${XOLOTL_UTIL_HEADER_DIR}/StepSequence.h
    ${XOLOTL_UTIL_HEADER_DIR}/LinearStepSequence.h
    ${XOLOTL_UTIL_HEADER_DIR}/GrowthFactorStepSequence.h
    ${XOLOTL_UTIL_HEADER_DIR}/Log.h
    ${XOLOTL_UTIL_HEADER_DIR}/MathUtils.h
    ${XOLOTL_UTIL_HEADER_DIR}/MPIUtils.h
    ${XOLOTL_UTIL_HEADER_DIR}/NotImplementedError.h
    ${XOLOTL_UTIL_HEADER_DIR}/StreamTokenizer.h
    ${XOLOTL_UTIL_HEADER_DIR}/TimeStepper.h
    ${XOLOTL_UTIL_HEADER_DIR}/Tokenizer.h
    ${XOLOTL_UTIL_HEADER_DIR}/TupleUtils.h
    ${XOLOTL_UTIL_HEADER_DIR}/RandomNumberGenerator.h
)

set(XOLOTL_UTIL_SOURCES
    ${XOLOTL_UTIL_SOURCE_DIR}/GrowthFactorStepSequence.cpp
    ${XOLOTL_UTIL_SOURCE_DIR}/LinearStepSequence.cpp
    ${XOLOTL_UTIL_SOURCE_DIR}/Log.cpp
    ${XOLOTL_UTIL_SOURCE_DIR}/MPIUtils.cpp
    ${XOLOTL_UTIL_SOURCE_DIR}/TimeStepper.cpp
)

add_library(xolotlUtil SHARED
    ${XOLOTL_UTIL_SOURCES}
    ${XOLOTL_UTIL_HEADERS}
)
target_link_libraries(xolotlUtil PUBLIC
    MPI::MPI_CXX
    Boost::log_setup
    Boost::log
    Boost::dynamic_linking
    Kokkos::kokkos
)
target_include_directories(xolotlUtil PUBLIC
    $<BUILD_INTERFACE:${XOLOTL_UTIL_INCLUDE_DIR}>
    $<BUILD_INTERFACE:${XOLOTL_BINARY_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(xolotlUtil PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS xolotlUtil EXPORT Xolotl LIBRARY DESTINATION lib)
