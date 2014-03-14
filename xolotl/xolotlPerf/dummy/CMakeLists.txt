#dummy

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp)



#install(FILES ${HEADERS} DESTINATION include)
