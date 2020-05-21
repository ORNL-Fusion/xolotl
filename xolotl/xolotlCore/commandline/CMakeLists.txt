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

#Add a library to hold the command line handling code
add_library(xolotlCL STATIC ${SRC})
target_link_libraries(xolotlCL PUBLIC
    xolotlIO
    xolotlPerf
    Boost::program_options
    MPI::MPI_CXX
)
target_include_directories(xolotlCL PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
