#Set the package name
SET(PACKAGE_NAME "xolotl.diffusion")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Diffusion")
#Set the library name
SET(LIBRARY_NAME "xolotlDiffusion")

#Collect all header filenames in this project
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Add a library
add_library(xolotlDiffusion STATIC ${SRC})
target_link_libraries(xolotlDiffusion PUBLIC
    xolotlAdvection
    xolotlReactants
)
target_include_directories(xolotlDiffusion PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
