list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/temperature/HeatEquationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/ITemperatureHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/ConstantHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/GradientHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/TemperatureHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/temperature/ProfileHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/temperature/ConstantHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/temperature/GradientHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/temperature/HeatEquationHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/temperature/ProfileHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/temperature/TemperatureHandler.cpp
)
