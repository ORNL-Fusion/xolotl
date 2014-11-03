#Set the package name
SET(PACKAGE_NAME "xolotl.factory")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Factory")
#Set the library name
SET(LIBRARY_NAME "xolotlFactory")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS material/*.h vizHandler/*.h temperatureHandler/*.h)

#Grab all of the source files
file(GLOB SRC material/*.cpp vizHandler/*.cpp temperatureHandler/*.cpp)

#Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/material)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/advection)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/flux)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/diffusion)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/vizHandler)
include_directories(${CMAKE_BINARY_DIR}/xolotlFactory)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/dummy)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/temperatureHandler)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/temperature)

# Find EAVL
FIND_PACKAGE(EAVL)
# Find MESA
FIND_PACKAGE(MESA)
# If they are present
if(EAVL_FOUND AND MESA_FOUND)
    include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/standard)
endif(EAVL_FOUND AND MESA_FOUND)

#Add a library to hold the factory code
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlViz)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
