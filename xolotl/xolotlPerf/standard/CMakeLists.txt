#Set the package name
SET(PACKAGE_NAME "xolotl.performance.standard")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp)

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)
include_directories(${GPTL_INCLUDE_DIRS})
include_directories(${PAPI_INCLUDE_DIRS})

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)