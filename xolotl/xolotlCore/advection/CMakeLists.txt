#Set the package name
SET(PACKAGE_NAME "xolotl.advection")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Advection")
#Set the library name
SET(LIBRARY_NAME "xolotlAdvection")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Add a library
add_library(xolotlAdvection STATIC ${SRC})
target_link_libraries(xolotlAdvection PUBLIC
    xolotlCL
    xolotlReactants
    xolotlPerf
)
target_include_directories(xolotlAdvection PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
