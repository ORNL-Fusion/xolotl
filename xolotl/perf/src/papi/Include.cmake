list(APPEND XOLOTL_PERF_HEADERS
    ${XOLOTL_PERF_HEADER_DIR}/papi/PAPIHandlerRegistry.h
    ${XOLOTL_PERF_HEADER_DIR}/papi/PAPIHardwareCounter.h
    ${XOLOTL_PERF_HEADER_DIR}/papi/PAPITimer.h
)

list(APPEND XOLOTL_PERF_SOURCES
    ${XOLOTL_PERF_SOURCE_DIR}/papi/PAPIHandlerRegistry.cpp
    ${XOLOTL_PERF_SOURCE_DIR}/papi/PAPIHardwareCounter.cpp
    ${XOLOTL_PERF_SOURCE_DIR}/papi/PAPITimer.cpp
)
