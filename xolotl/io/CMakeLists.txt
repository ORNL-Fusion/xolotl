## NOTE: internal library (?)

set(XOLOTL_IO_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_IO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_IO_HEADER_DIR ${XOLOTL_IO_INCLUDE_DIR}/xolotl/io)

set(XOLOTL_IO_HEADERS
    ${XOLOTL_IO_HEADER_DIR}/HDF5Exception.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5File.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileAttribute.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileDataSet.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileDataSpace.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileGroup.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5FileType.h
    ${XOLOTL_IO_HEADER_DIR}/HDF5Object.h
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
    ${XOLOTL_IO_SOURCE_DIR}/XFile.cpp
)

add_library(xolotlIO SHARED
    ${XOLOTL_IO_SOURCES}
    ${XOLOTL_IO_HEADERS}
)
target_link_libraries(xolotlIO PUBLIC
    xolotlCore
    xolotlUtil
    ${HDF5_LIBRARIES}
)
target_include_directories(xolotlIO PUBLIC
    ${HDF5_INCLUDE_DIR}
    $<BUILD_INTERFACE:${XOLOTL_IO_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:include>
)
set_target_properties(xolotlIO PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(TARGETS xolotlIO EXPORT Xolotl LIBRARY DESTINATION lib)
