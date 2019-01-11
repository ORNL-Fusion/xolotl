# use VTKM_DIR to setup the options that cmake's find VTKm needs
set(VTKm_DIR ${VTKM_DIR}/lib/cmake/vtkm-1.3)
MESSAGE(STATUS "Looking for VTKm using VTKM_DIR = ${VTKm_DIR}")

find_package(VTKm REQUIRED QUIET)
set(VTKm_LIBRARIES vtkm_cont vtkm_rendering)
set(VTKm_INCLUDE_DIRS "${VTKm_DIR}/include/vtkm-1.3 ${VTKm_DIR}/vtkm/thirdparty/taotuple")
set(VTKM_FOUND TRUE)

message(STATUS "VTKm includes = ${VTKm_INCLUDE_DIRS}")
message(STATUS "VTKm libraries = ${VTKm_LIBRARIES}")
