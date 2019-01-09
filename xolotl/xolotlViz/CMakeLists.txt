#Set the package name
SET(PACKAGE_NAME "xolotl.viz")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Visualization")
#Set the library name
SET(LIBRARY_NAME "xolotlViz")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dataprovider)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/labelprovider)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dummy)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Always include the Dummy (stub) classes.
file(GLOB HEADERS *.h dummy/*.h dataprovider/*.h labelprovider/*.h)
file(GLOB SRC *.cpp dummy/*.cpp dataprovider/*.cpp)

# Find VTKm
if(VTKM_DIR)
    include(${CMAKE_SOURCE_DIR}/cmake-modules/SetupVTKm.cmake)
else()
    message(STATUS "Visualization support needs explicit VTKm_DIR.")
endif()

if(VTKM_FOUND)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/standard)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/standard/plot)
    include_directories(${VTKm_INCLUDE_DIRS})

    # Add the standard classes to the HEADERS and SRCS.
    file(GLOB STDHEADERS standard/*.h standard/plot/*.h)
    set(HEADERS ${HEADERS} ${STDHEADERS})
    file(GLOB STDSRC standard/*.cpp standard/plot/*.cpp)
    set(SRC ${SRC} ${STDSRC})

    set(MAYBE_VTKM ${VTKm_LIBRARIES})

    # Ensure that our config file indicates we have built the
    # standard visualization classes in our visualization infrastructure.
    set(HAVE_VIZLIB_STD 1)
endif(VTKM_FOUND)


#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
message('${SRC}')
message('${HEADERS}')
target_link_libraries(${LIBRARY_NAME} ${MAYBE_VTKM})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

# Configure file for the performance library.
# Note: must do this after all the other checks, or else the
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_SOURCE_DIR}/xolotlFactory/vizHandler/XolotlConfigViz.h.in" "${CMAKE_BINARY_DIR}/xolotlFactory/XolotlConfigViz.h" )
