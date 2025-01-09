find_package(Git REQUIRED)

message(STATUS "Setup PETSc")

## Set directories
set(__external_src_dir ${CMAKE_SOURCE_DIR}/external)
set(__external_bin_dir ${CMAKE_BINARY_DIR}/external)

set(__petsc_src_dir ${__external_src_dir}/petsc)
set(__petsc_bin_dir ${__external_bin_dir}/petsc_build)

execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${__petsc_bin_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${__petsc_bin_dir})

if(NOT EXISTS ${__petsc_src_dir}/configure)
    message(STATUS "    checkout")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init external/petsc
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_FILE "${__external_bin_dir}/petsc_clone.out"
        ERROR_FILE "${__external_bin_dir}/petsc_clone.out"
    )
endif()

## Check for dependencies
#### HDF5
option(Xolotl_BUILD_HDF5 "Have the PETSc build system build HDF5" OFF)
if(NOT Xolotl_BUILD_HDF5)
    set(HDF5_PREFER_PARALLEL ON)
    find_package(HDF5 QUIET)
endif()
#### Boost
option(Xolotl_BUILD_BOOST "Have the PETSc build system build Boost" OFF)
if(NOT Xolotl_BUILD_BOOST)
    find_package(Boost QUIET COMPONENTS
        log_setup
        log
        program_options
    )
endif()
#### BLAS/LAPACK
option(Xolotl_BUILD_LAPACK "Have the PETSc build system build LAPACK" OFF)
if(NOT Xolotl_BUILD_LAPACK)
    find_package(LAPACK QUIET)
endif()

## Set options
set(__script_dir ${CMAKE_SOURCE_DIR}/scripts)
set(__build_opts
    --skip-pull
    --prefix=${__external_bin_dir}/petsc_install
)

option(Xolotl_BUILD_PETSC_DEBUG
    "Enable debugging symbols for petsc, kokkos, etc."
    OFF
)
if(Xolotl_BUILD_PETSC_DEBUG)
    list(APPEND __build_opts --debug)
    message(STATUS "    - enable debugging")
endif()

option(Xolotl_ENABLE_CUDA "Enable CUDA backend for kokkos, etc." OFF)
option(Xolotl_ENABLE_OPENMP "Enable OpenMP backend for kokkos" OFF)
if(Xolotl_ENABLE_CUDA)
    list(APPEND __build_opts --cuda)
    message(STATUS "    - using CUDA backend")
    if(Xolotl_CUDA_SM)
        list(APPEND __build_opts --cuda-sm=${Xolotl_CUDA_SM})
    else()
        execute_process(
            COMMAND nvidia-smi
            OUTPUT_VARIABLE __smi_out
            ERROR_VARIABLE __smi_out
            RESULT_VARIABLE __smi_avail
        )
        if(NOT ${__smi_avail} EQUAL 0)
            message(FATAL_ERROR "
                Unable to determine CUDA SM version (compute capability).
                Please provide this in Xolotl_CUDA_SM.
                "
            )
        endif()
    endif()
elseif(Xolotl_ENABLE_OPENMP)
    list(APPEND __build_opts --openmp)
    message(STATUS "    - using OpenMP backend")
else()
    message(STATUS "    - using Serial backend")
endif()

if(Xolotl_KOKKOS_VERSION)
    list(APPEND __build_opts --kokkos-version=${Xolotl_KOKKOS_VERSION})
endif()

## Include dependencies if necessary
#### HDF5
if(NOT HDF5_FOUND)
    set(Xolotl_BUILD_HDF5 ON)
endif()
if(Xolotl_BUILD_HDF5)
    list(APPEND __build_opts --get-hdf5)
    message(STATUS "    - build HDF5")
    set(HDF5_ROOT "${__external_bin_dir}/petsc_install")
endif()
#### Boost
if(NOT Boost_FOUND)
    set(Xolotl_BUILD_BOOST ON)
endif()
if(Xolotl_BUILD_BOOST)
    list(APPEND __build_opts --get-boost)
    message(STATUS "    - build Boost")
    if(NOT Xolotl_BUILD_PETSC_DEBUG)
        set(Boost_USE_DEBUG_RUNTIME OFF CACHE INTERNAL "")
    endif()
endif()
#### LAPACK
if(NOT LAPACK_FOUND)
    set(Xolotl_BUILD_LAPACK ON)
endif()
if(Xolotl_BUILD_LAPACK)
    list(APPEND __build_opts --get-lapack)
    message(STATUS "    - build BLAS/LAPACK")
endif()

## Perform build
set(__output_file "${__external_bin_dir}/petsc_build.out")
message(STATUS "    build (for output, follow ${__output_file})")
set(__command bash ${__script_dir}/build_petsc.sh ${__build_opts})
if(Xolotl_BUILD_PETSC_DRY_RUN)
    string(REPLACE ";" " " __command_str "${__command}")
    message(STATUS "Script Command:")
    message(STATUS "${__command_str}")
    execute_process(
        COMMAND ${__command} --dry-run
        WORKING_DIRECTORY "${__petsc_src_dir}"
    )
    message(FATAL_ERROR "Exiting")
endif()
execute_process(
    COMMAND ${__command}
    WORKING_DIRECTORY "${__petsc_src_dir}"
    OUTPUT_FILE "${__output_file}"
    ERROR_FILE "${__output_file}"
    RESULT_VARIABLE __build_ret
)
if(NOT ${__build_ret} EQUAL 0)
    message(FATAL_ERROR "
        Failed to build PETSc
        See \"${__output_file}\"
        "
    )
endif()

list(APPEND CMAKE_PREFIX_PATH "${__external_bin_dir}/petsc_install")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} CACHE PATH "" FORCE)

## Don't build when re-running CMake unless the user specifies this again
set(Xolotl_BUILD_PETSC OFF CACHE PATH "" FORCE)
