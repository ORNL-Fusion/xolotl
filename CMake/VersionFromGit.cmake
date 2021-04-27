include(CMakeParseArguments)

function(version_from_git)
    # Parse arguments
    set(options OPTIONAL FAST)
    set(oneValueArgs
        GIT_EXECUTABLE
        LOG
        TIMESTAMP
    )
    set(multiValueArgs)
    cmake_parse_arguments(ARG
        "${options}"
        "${oneValueArgs}"
        "${multiValueArgs}"
        ${ARGN}
    )

    if(DEFINED ARG_GIT_EXECUTABLE)
        set( GIT_EXECUTABLE "${ARG_GIT_EXECUTABLE}")
    else()
        # Find Git or bail out
        find_package(Git)
        if(NOT GIT_FOUND)
            message(FATAL_ERROR "Git not found")
        endif(NOT GIT_FOUND)
    endif()

    # Git describe
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --tags
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE git_result
        OUTPUT_VARIABLE git_describe
        ERROR_VARIABLE git_error
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
    )
    if(NOT git_result EQUAL 0)
        message(FATAL_ERROR "${git_error}")
    endif()

    # Get Git tag
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --tags --abbrev=0
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE git_result
        OUTPUT_VARIABLE git_tag
        ERROR_VARIABLE git_error
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_STRIP_TRAILING_WHITESPACE
    )
    if( NOT git_result EQUAL 0 )
        message(FATAL_ERROR "${git_error}")
    endif()

    if(git_tag MATCHES "^v([0-9]+)[.]([0-9]+)[.]([0-9]+)[+]?([.0-9A-Za-z-]+)?$")
        set(version_major "${CMAKE_MATCH_1}")
        set(version_minor "${CMAKE_MATCH_2}")
        set(version_patch "${CMAKE_MATCH_3}")
        set(version_label "${CMAKE_MATCH_4}")
    else()
        message(FATAL_ERROR "Git tag isn't valid version: [${git_tag}]")
    endif()

    if("${git_tag}" STREQUAL "${git_describe}")
        set(git_at_a_tag ON)
    endif()

    if(NOT git_at_a_tag)
        string(REPLACE "${git_tag}" "" post_tag_desc "${git_describe}")
        # Extract the Git revision and hash
        if(post_tag_desc MATCHES "^[-]([0-9]+)[-][g]([0-9a-f]+)$")
            set(revision "${CMAKE_MATCH_1}")
            set(git_hash "${CMAKE_MATCH_2}")
        endif()
    endif()

    # Construct the version variables
    set(version ${version_major}.${version_minor}.${version_patch})
    set(semver ${version})

    if(version_label MATCHES ".+")
        set(semver "${semver}+${version_label}")
    endif()

    # Log the results
    if(ARG_LOG)
        message(STATUS "Git Version Info:
            Git describe:   [${git_describe}]
            Git tag:        [${git_tag}]
            Version:        [${version}]
            Label:          [${version_label}]
            Revision:       [${revision}]
            Git hash:       [${git_hash}]
            SemVer:         [${semver}]"
        )
    endif(ARG_LOG)

    # Set parent scope variables
    set(GIT_TAG ${git_tag} PARENT_SCOPE)
    set(SEMVER ${semver} PARENT_SCOPE)
    set(VERSION ${version} PARENT_SCOPE)
    set(VERSION_MAJOR ${version_major} PARENT_SCOPE)
    set(VERSION_MINOR ${version_minor} PARENT_SCOPE)
    set(VERSION_PATCH ${version_patch} PARENT_SCOPE)

endfunction(version_from_git)

version_from_git(LOG ON)
