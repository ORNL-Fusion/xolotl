#Add the core subdirectories
add_subdirectory(commandline)
add_subdirectory(io)
add_subdirectory(reactants)
add_subdirectory(diffusion)
add_subdirectory(advection)
add_subdirectory(flux)
add_subdirectory(temperature)
add_subdirectory(modifiedreaction)

# Copy the headers for install(
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION ${PROJECT_BINARY_DIR}/include)