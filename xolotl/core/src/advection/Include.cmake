list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/advection/AdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/DummyAdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/IAdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/SurfaceAdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/TungstenAdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/W100AdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/W110AdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/W111AdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/W211AdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/XGBAdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/YGBAdvectionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/advection/ZGBAdvectionHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/advection/SurfaceAdvectionHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/advection/TungstenAdvectionHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/advection/XGBAdvectionHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/advection/YGBAdvectionHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/advection/ZGBAdvectionHandler.cpp
)
