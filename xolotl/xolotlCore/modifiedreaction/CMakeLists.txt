#Set the package name
SET(PACKAGE_NAME "xolotl.modified")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Modified Reaction")
#Set the library name
SET(LIBRARY_NAME "xolotlModified")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS trapmutation/*.h heterogeneousnucleation/*.h)

#Grab all of the source files
file(GLOB SRC trapmutation/*.cpp heterogeneousnucleation/*.cpp)

#Add a library
add_library(xolotlModified STATIC ${SRC})
target_link_libraries(xolotlModified PUBLIC
    xolotlAdvection
    xolotlCL
    xolotlReactants
    xolotlPerf
)
target_include_directories(xolotlModified PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/heterogeneousnucleation
    ${CMAKE_CURRENT_SOURCE_DIR}/trapmutation
)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
