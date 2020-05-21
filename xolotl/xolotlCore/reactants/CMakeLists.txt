#Set the package name
SET(PACKAGE_NAME "xolotl.reactants")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Reactants")
#Set the library name
SET(LIBRARY_NAME "xolotlReactants")

#Collect all header filenames in this project
#and glob them in HEADERS
file(GLOB HEADERS *.h psiclusters/*.h neclusters/*.h feclusters/*.h alloyclusters/*.h)

#Collect all of the cpp files in this folder
#and the psiclusters folder and glob them in SRC
file(GLOB SRC *.cpp psiclusters/*.cpp neclusters/*.cpp feclusters/*.cpp alloyclusters/*.cpp)

#Add the library
add_library(xolotlReactants STATIC ${SRC})
target_link_libraries(xolotlReactants PUBLIC
    xolotlPerf
    xolotlIO
    xolotlCL
    plsm
    MPI::MPI_CXX
)
target_include_directories(xolotlReactants PUBLIC
    ${CMAKE_SOURCE_DIR}/xolotlCore
    ${CMAKE_SOURCE_DIR}/xolotlCore/reactants
    ${CMAKE_CURRENT_SOURCE_DIR}/psiclusters
    ${CMAKE_CURRENT_SOURCE_DIR}/neclusters
    ${CMAKE_CURRENT_SOURCE_DIR}/feclusters
    ${CMAKE_CURRENT_SOURCE_DIR}/alloyclusters
)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS xolotlReactants DESTINATION lib)
