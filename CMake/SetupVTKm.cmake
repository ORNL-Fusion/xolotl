# use VTKM_DIR to setup the options that cmake's find VTKm needs
if(NOT VTKM_DIR)
    if(VTKm_DIR)
        set(VTKM_DIR ${VTKm_DIR})
    else()
        message(FATAL_ERROR "VTKm support needs explicit VTKM_DIR")
    endif()
endif()

message(STATUS "Looking for VTKm using VTKM_DIR = ${VTKM_DIR}")

# use VTKM_DIR to setup the options that cmake's find VTKm needs
file(GLOB VTKm_DIR "${VTKM_DIR}/lib/cmake/vtkm-*")
if(NOT VTKm_DIR)
    message(FATAL_ERROR "Failed to find VTKm at VTKM_DIR=${VTKM_DIR}/lib/cmake/vtk-*")
endif()

find_package(VTKm REQUIRED)
set(VTKm_LIBRARIES vtkm_worklet vtkm_cont vtkm_rendering)
set(VTKM_FOUND TRUE)

message(STATUS "VTKm includes = ${VTKm_INCLUDE_DIRS}")
message(STATUS "VTKm libraries = ${VTKm_LIBRARIES}")
