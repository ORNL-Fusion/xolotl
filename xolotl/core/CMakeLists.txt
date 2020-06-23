set(XOLOTL_CORE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(XOLOTL_CORE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XOLOTL_CORE_INCLUDE_DIR ${XOLOTL_CORE_INCLUDE_DIR} PARENT_SCOPE)
set(XOLOTL_CORE_HEADER_DIR ${XOLOTL_CORE_INCLUDE_DIR}/xolotl/core)

set(XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/Constants.h
)

include(src/network/Include.cmake)
include(src/diffusion/Include.cmake)
include(src/advection/Include.cmake)
include(src/flux/Include.cmake)
include(src/temperature/Include.cmake)
include(src/modified/Include.cmake)

add_library(xolotlCore STATIC
    ${XOLOTL_CORE_SOURCES}
    ${XOLOTL_CORE_HEADERS}
    ${XOLOTL_CORE_TPP_FILES}
)
target_link_libraries(xolotlCore PUBLIC
    xolotlPerf
    xolotlOptions
    xolotlUtil
    plsm
    MPI::MPI_CXX
)
target_include_directories(xolotlCore PUBLIC ${XOLOTL_CORE_INCLUDE_DIR})

install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/advection/IAdvectionHandler.h
    DESTINATION include/xolotl/core/advection
)
install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/IDiffusionHandler.h
    DESTINATION include/xolotl/core/diffusion
)
install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/flux/IFluxHandler.h
    DESTINATION include/xolotl/core/flux
)
install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/modified/ITrapMutationHandler.h
    DESTINATION include/xolotl/core/modified
)
install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/temperature/ITemperatureHandler.h
    DESTINATION include/xolotl/core/temperature
)
install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/network/Cluster.h
    ${XOLOTL_CORE_HEADER_DIR}/network/ClusterData.h
    ${XOLOTL_CORE_HEADER_DIR}/network/IReactionNetwork.h
    ${XOLOTL_CORE_HEADER_DIR}/network/ReactionNetworkTraits.h
    ${XOLOTL_CORE_HEADER_DIR}/network/SpeciesEnumSequence.h
    ${XOLOTL_CORE_HEADER_DIR}/network/EnumSequence.h
    DESTINATION include/xolotl/core/network
)
install(FILES
    ${XOLOTL_CORE_HEADER_DIR}/network/detail/MemorySpace.h
    ${XOLOTL_CORE_HEADER_DIR}/network/detail/ReactionData.h
    DESTINATION include/xolotl/core/network/detail
)
