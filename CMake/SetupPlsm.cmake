option(Xolotl_USE_PLSM_DEVELOP "" FALSE)
if(EXISTS ${plsm_DIR})
    if(DEFINED CACHE{XOLOTL_USE_PLSM_DEVELOP})
        if(Xolotl_USE_PLSM_DEVELOP AND XOLOTL_USE_PLSM_DEVELOP)
            set(__return TRUE)
        elseif(NOT Xolotl_USE_PLSM_DEVELOP AND NOT XOLOTL_USE_PLSM_DEVELOP)
            set(__return TRUE)
        else()
            set(__return FALSE)
            set(__do_src_update TRUE)
        endif()
    endif()
    if(DEFINED CACHE{XOLOTL_USE_64BIT_INDEX_TYPE})
        if(Xolotl_USE_64BIT_INDEX_TYPE AND XOLOTL_USE_64BIT_INDEX_TYPE)
        elseif(NOT Xolotl_USE_64BIT_INDEX_TYPE AND NOT
                XOLOTL_USE_64BIT_INDEX_TYPE)
        else()
            set(__return FALSE)
        endif()
    endif()
    if(__return)
        return()
    endif()
endif()

find_package(Git REQUIRED)

message(STATUS "Setup plsm")

set(__external_src_dir ${CMAKE_SOURCE_DIR}/external)
set(__external_bin_dir ${CMAKE_BINARY_DIR}/external)
set(__plsm_src_dir ${__external_src_dir}/plsm)
set(__plsm_bin_dir ${__external_bin_dir}/plsm_build)

execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${__plsm_bin_dir})
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${__plsm_bin_dir})

if(NOT EXISTS ${__plsm_src_dir})
    message(STATUS "    checkout")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_FILE "${__external_bin_dir}/plsm_clone.out"
        ERROR_FILE "${__external_bin_dir}/plsm_clone.out"
    )
    execute_process(
        COMMAND ${GIT_EXECUTABLE} fetch --all
        WORKING_DIRECTORY "${__plsm_src_dir}"
        OUTPUT_FILE "${__external_bin_dir}/plsm_clone.out"
        ERROR_FILE "${__external_bin_dir}/plsm_clone.out"
    )
endif()

# Should we update the source to a different version
set(XOLOTL_USE_PLSM_DEVELOP ${Xolotl_USE_PLSM_DEVELOP} CACHE INTERNAL "")
if(__do_src_update)
    message(STATUS "    checkout")
    if(${Xolotl_USE_PLSM_DEVELOP})
        execute_process(
            COMMAND ${GIT_EXECUTABLE} checkout origin/develop
            WORKING_DIRECTORY "${__plsm_src_dir}"
            OUTPUT_FILE "${__external_bin_dir}/plsm_clone.out"
            ERROR_FILE "${__external_bin_dir}/plsm_clone.out"
        )
    else()
        # Get original submodule commit and checkout
        execute_process(
            COMMAND ${GIT_EXECUTABLE} ls-tree HEAD:external
            WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
            OUTPUT_VARIABLE __external_ls_output
        )
        separate_arguments(__external_ls_output
            NATIVE_COMMAND
            ${__external_ls_output}
        )
        list(GET __external_ls_output 2 __ref_commit)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} checkout ${__ref_commit}
            WORKING_DIRECTORY "${__plsm_src_dir}"
            OUTPUT_FILE "${__external_bin_dir}/plsm_clone.out"
            ERROR_FILE "${__external_bin_dir}/plsm_clone.out"
        )
    endif()
endif()

set(__plsm_opts
    -DKokkos_DIR=${Kokkos_DIR}
    -DCMAKE_INSTALL_PREFIX=${__plsm_bin_dir}/install
    -DBUILD_TESTING=OFF
    -DPLSM_USE_64BIT_INDEX_TYPE=${Xolotl_USE_64BIT_INDEX_TYPE}
)
message(STATUS "    configure")
set(__output_file "${__external_bin_dir}/plsm_configure.out")
execute_process(
    COMMAND ${CMAKE_COMMAND} ${__plsm_opts} ${__plsm_src_dir}
    WORKING_DIRECTORY "${__plsm_bin_dir}"
    OUTPUT_FILE "${__output_file}"
    ERROR_FILE "${__output_file}"
    RESULT_VARIABLE __config_ret
)
if(NOT ${__config_ret} EQUAL 0)
    message(FATAL_ERROR "
        Failed to configure plsm
        See \"${__output_file}\"
        "
    )
endif()

message(STATUS "    install")
set(__output_file "${__external_bin_dir}/plsm_build.out")
execute_process(
    COMMAND ${CMAKE_COMMAND} --build . --target install
    WORKING_DIRECTORY "${__plsm_bin_dir}"
    OUTPUT_FILE "${__output_file}"
    ERROR_FILE "${__output_file}"
    RESULT_VARIABLE __build_ret
)
if(NOT ${__build_ret} EQUAL 0)
    message(FATAL_ERROR "
        Failed to install plsm
        See \"${__output_file}\"
        "
    )
endif()

set(plsm_DIR ${__plsm_bin_dir}/install)
