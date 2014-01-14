#Set the package name
SET(PACKAGE_NAME "xolotl.io")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl I/O")
#Set the library name
SET(LIBRARY_NAME "xolotlIO")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Include the MPI directories
include_directories(${MPI_INCLUDE_PATH})

#Add a library to hold the IO code (mostly MPI)
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} ${MPI_LIBRARIES})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
