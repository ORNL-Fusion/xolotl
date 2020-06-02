# Find VTKm
if(VTKM_DIR)
    include(SetupVTKm)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

if(VTKM_FOUND)
    # Ensure that our config file indicates we have built the
    # standard visualization classes in our visualization infrastructure.
    set(HAVE_VIZLIB_STD 1)
endif(VTKM_FOUND)
execute_process(COMMAND
    ${CMAKE_COMMAND} -E make_directory ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/viz
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
    "${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/viz/config.h"
)

set(XOLOTL_VIZ_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_VIZ_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_VIZ_HEADER_DIR ${XOLOTL_VIZ_INCLUDE_DIR}/xolotl/viz)

set(XOLOTL_VIZ_HEADERS
    ${XOLOTL_BINARY_INCLUDE_DIR}/xolotl/viz/config.h
    ${XOLOTL_VIZ_HEADER_DIR}/IPlot.h
    ${XOLOTL_VIZ_HEADER_DIR}/IVizHandlerRegistry.h
    ${XOLOTL_VIZ_HEADER_DIR}/LabelProvider.h
    ${XOLOTL_VIZ_HEADER_DIR}/PlottingStyle.h
    ${XOLOTL_VIZ_HEADER_DIR}/PlotType.h
    ${XOLOTL_VIZ_HEADER_DIR}/Point.h
)

include(src/dataprovider/Include.cmake)
include(src/dummy/Include.cmake)
if(VTKM_FOUND)
    include(src/standard/Include.cmake)
endif(VTKM_FOUND)

add_library(xolotlViz STATIC ${XOLOTL_VIZ_SOURCES} ${XOLOTL_VIZ_HEADERS})
target_link_libraries(xolotlViz PUBLIC xolotlCore)
target_include_directories(xolotlViz PUBLIC
    ${XOLOTL_VIZ_INCLUDE_DIR}
    ${XOLOTL_BINARY_INCLUDE_DIR}
)
if(VTKM_FOUND)
    target_include_directories(xolotlViz PUBLIC ${VTKm_INCLUDE_DIRS})
    target_link_libraries(xolotlViz PUBLIC ${VTKm_LIBRARIES})
endif(VTKM_FOUND)

install(FILES ${XOLOTL_VIZ_HEADERS} DESTINATION include)
install(TARGETS xolotlViz DESTINATION lib)