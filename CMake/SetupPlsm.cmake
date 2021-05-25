option(Xolotl_USE_PLSM_DEVELOP "" FALSE)
if(EXISTS ${plsm_DIR})
    if(DEFINED CACHE{XOLOTL_USE_PLSM_DEVELOP})
        if(Xolotl_USE_PLSM_DEVELOP AND XOLOTL_USE_PLSM_DEVELOP)
            return()
        elseif(NOT Xolotl_USE_PLSM_DEVELOP AND NOT XOLOTL_USE_PLSM_DEVELOP)
            return()
        else()
            set(__do_src_update TRUE)
        endif()
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

function(init_plsm_submodule)
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
endfunction()

if(NOT EXISTS ${__plsm_src_dir})
    message(STATUS "    checkout")
    init_plsm_submodule()
endif()

# Should we update the source to a different version
set(XOLOTL_USE_PLSM_DEVELOP ${Xolotl_USE_PLSM_DEVELOP} CACHE INTERNAL "")
if(__do_src_update)
    if(${Xolotl_USE_PLSM_DEVELOP})
        execute_process(
            COMMAND ${GIT_EXECUTABLE} checkout origin/develop
            WORKING_DIRECTORY "${__plsm_src_dir}"
            OUTPUT_FILE "${__external_bin_dir}/plsm_clone.out"
            ERROR_FILE "${__external_bin_dir}/plsm_clone.out"
        )
    else()
        init_plsm_submodule()
    endif()
endif()

set(__plsm_opts
    -DKokkos_DIR=${Kokkos_DIR}
    -DCMAKE_INSTALL_PREFIX=${__plsm_bin_dir}/install
    -DBUILD_TESTING=OFF
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
