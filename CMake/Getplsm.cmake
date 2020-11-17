set(PLSM_URL https://github.com/ORNL-Fusion/plsm.git)

find_package(Git REQUIRED)

message(STATUS "Get plsm")

set(__external_dir ${CMAKE_BINARY_DIR}/external)
set(__plsm_src_dir ${__external_dir}/plsm_src)
set(__plsm_bin_dir ${__external_dir}/plsm)
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${__plsm_bin_dir})
if(NOT EXISTS ${__plsm_src_dir})
    message(STATUS "    checkout")
    execute_process(COMMAND ${GIT_EXECUTABLE} clone ${PLSM_URL} plsm_src
        WORKING_DIRECTORY ${__external_dir}
        OUTPUT_FILE ${__external_dir}/plsm_clone.out
        ERROR_FILE ${__external_dir}/plsm_clone.out
    )
endif()

set(__plsm_opts
    -DKokkos_DIR=${Kokkos_DIR}
    -DCMAKE_INSTALL_PREFIX=${__plsm_bin_dir}/install
    -DBUILD_TESTING=OFF
)
message(STATUS "    configure")
execute_process(COMMAND ${CMAKE_COMMAND} ${__plsm_opts} ${__plsm_src_dir}
    WORKING_DIRECTORY ${__plsm_bin_dir}
    OUTPUT_FILE ${__external_dir}/plsm_configure.out
    ERROR_FILE ${__external_dir}/plsm_configure.out
)

message(STATUS "    install")
execute_process(COMMAND ${CMAKE_COMMAND} --build . --target install
    WORKING_DIRECTORY ${__plsm_bin_dir}
    OUTPUT_FILE ${__external_dir}/plsm_install.out
    ERROR_FILE ${__external_dir}/plsm_install.out
)

set(plsm_DIR ${__plsm_bin_dir}/install)
