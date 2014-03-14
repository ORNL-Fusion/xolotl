#standard

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)

# Find GPTL
FIND_PACKAGE(GPTL)

# Find PAPI
FIND_PACKAGE(PAPI)

if(GPTL_FOUND AND PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")
    message(STATUS "GPTL includes = ${GPTL_INCLUDE_DIRS}")
    message(STATUS "GPTL libraries = ${GPTL_LIBRARIES}")
    
    #Collect all header filenames in this project 
    #and glob them in HEADERS
    file(GLOB HEADERS *.h)
    
    message(STATUS "HEADERS = ${HEADERS}")

    #Collect all of the cpp files in this folder 
    #and glob them in SRC
    file(GLOB SRC *.cpp)

    #Add the library
    add_library(${LIBRARY_NAME} SHARED ${SRC})
    target_link_libraries(${LIBRARY_NAME} xolotlPerf ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})
    
endif(GPTL_FOUND AND PAPI_FOUND)

install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)