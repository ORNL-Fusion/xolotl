#Set the package name
SET(PACKAGE_NAME "xolotl.commandline")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Command Line Option Handling")
#Set the library name
SET(LIBRARY_NAME "xolotlCL")

#Collect all header filenames in this project
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

# Where to look for include files
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/io)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/)

#Add a library to hold the command line handling code
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME} ${Boost_LIBRARIES} MPI::MPI_CXX)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
