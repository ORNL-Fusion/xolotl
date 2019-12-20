list(APPEND CMAKE_PREFIX_PATH ${plsm_DIR})
find_package(plsm REQUIRED)

#Add the core subdirectories
add_subdirectory(commandline)
add_subdirectory(io)
add_subdirectory(reactants)
add_subdirectory(diffusion)
add_subdirectory(advection)
add_subdirectory(flux)
add_subdirectory(temperature)
add_subdirectory(modifiedreaction)
