find_package(Git REQUIRED)

message(STATUS "Setup PETSc")

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

# Check for dependencies
set(HDF5_PREFER_PARALLEL ON)
find_package(HDF5 QUIET)
find_package(Boost QUIET COMPONENTS
    log_setup
    log
    program_options
)
find_package(LAPACK QUIET)

# Set options
set(__script_dir ${CMAKE_SOURCE_DIR}/scripts)
option(Xolotl_ENABLE_CUDA "Enable CUDA backend for kokkos, etc." OFF)
set(__build_opts
    --skip-pull
    --prefix=${__external_bin_dir}/petsc_install
)

if(CMAKE_BUILD_TYPE MATCHES "Debug")
    list(APPEND __build_opts --debug)
    message(STATUS "    - enable debugging")
endif()
if(Xolotl_ENABLE_CUDA)
    list(APPEND __build_opts --cuda)
    message(STATUS "    - use CUDA")
endif()
if(NOT HDF5_FOUND)
    list(APPEND __build_opts --get-hdf5)
    message(STATUS "    - build HDF5")
endif()
if(NOT Boost_FOUND)
    list(APPEND __build_opts --get-boost)
    message(STATUS "    - build Boost")
endif()
if(NOT LAPACK_FOUND)
    list(APPEND __build_opts --get-lapack)
    message(STATUS "    - build BLAS/LAPACK")
endif()
set(__output_file "${__external_bin_dir}/petsc_build.out")
message(STATUS "    build (for output, follow ${__output_file})")
execute_process(
    COMMAND bash ${__script_dir}/build_petsc.sh ${__build_opts}
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

list(APPEND CMAKE_PREFIX_PATH ${__external_bin_dir}/petsc_install)
