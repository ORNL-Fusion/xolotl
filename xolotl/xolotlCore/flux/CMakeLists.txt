#Set the package name
SET(PACKAGE_NAME "xolotl.flux")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Flux")
#Set the library name
SET(LIBRARY_NAME "xolotlFlux")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Add the libraries
add_library(xolotlFlux STATIC ${SRC})
target_link_libraries(xolotlFlux PUBLIC
    xolotlReactants
    xolotlCL
    xolotlPerf
)
target_include_directories(xolotlFlux PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
