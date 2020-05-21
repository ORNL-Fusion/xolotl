#Set the package name
SET(PACKAGE_NAME "xolotl.temperature")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Temperature")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

add_library(xolotlTemperature INTERFACE)
target_include_directories(xolotlTemperature INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
