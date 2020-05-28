## NOTE: internal library (?)

execute_process(COMMAND
    ${CMAKE_COMMAND} -E make_directory ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/io
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
    "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/io/config.h"
)

set(XOLOTL_IO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_IO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_IO_HEADER_DIR ${XOLOTL_IO_INCLUDE_DIR}/xolotl/io)

set(XOLOTL_IO_HEADERS
    ${XOLOTL_IO_HEADER_DIR}/Filesystem.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5Exception.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5File.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileAttribute.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileDataSet.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileDataSpace.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileGroup.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileType.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5Object.h
    ${XOLOTL_IO_HEADER_DIR}/MPIUtils.h
    ${XOLOTL_IO_HEADER_DIR}/TokenizedLineReader.h
    ${XOLOTL_IO_HEADER_DIR}/XFile.h
    ${XOLOTL_IO_HEADER_DIR}/XFileType.h
)

set(XOLOTL_IO_SOURCES
    ${XOLOTL_IO_SOURCE_DIR}/HDF5File.cpp
    ${XOLOTL_IO_SOURCE_DIR}/HDF5FileAttribute.cpp
    ${XOLOTL_IO_SOURCE_DIR}/HDF5FileDataSet.cpp
    ${XOLOTL_IO_SOURCE_DIR}/HDF5FileDataSpace.cpp
    ${XOLOTL_IO_SOURCE_DIR}/HDF5FileGroup.cpp
    ${XOLOTL_IO_SOURCE_DIR}/HDF5FileType.cpp
    ${XOLOTL_IO_SOURCE_DIR}/MPIUtils.cpp
    ${XOLOTL_IO_SOURCE_DIR}/XFile.cpp
)

add_library(xolotlIO STATIC ${XOLOTL_IO_SOURCES} ${XOLOTL_IO_HEADERS})
target_link_libraries(xolotlIO PUBLIC
    xolotlCore
    ${HDF5_LIBRARIES}
    MPI::MPI_CXX
    Kokkos::kokkos
    plsm
)
target_include_directories(xolotlIO PUBLIC
    ${XOLOTL_IO_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
)
if(Boost_FILESYSTEM_FOUND)
    target_link_libraries(xolotlIO PUBLIC Boost::filesystem)
endif(Boost_FILESYSTEM_FOUND)
