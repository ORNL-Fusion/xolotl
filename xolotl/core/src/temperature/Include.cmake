list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/temperature/HeatEquation1DHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/HeatEquation2DHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/HeatEquation3DHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/ITemperatureHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/TemperatureConstantHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/TemperatureGradientHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/TemperatureHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/TemperatureProfileHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/temperature/TemperatureHandler.cpp
)
