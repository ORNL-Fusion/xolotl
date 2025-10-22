list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/Diffusion1DHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/Diffusion2DHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/Diffusion3DHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/DiffusionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/DummyDiffusionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/IDiffusionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/diffusion/Spherical1DHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/diffusion/DiffusionHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/diffusion/Diffusion1DHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/diffusion/Diffusion2DHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/diffusion/Diffusion3DHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/diffusion/Spherical1DHandler.cpp
)
