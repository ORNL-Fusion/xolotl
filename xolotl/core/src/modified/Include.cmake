list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/modified/DummyTrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/ITrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/Sigma3TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/W100TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/W110TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/W111TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/W211TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/ISoretDiffusionHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modified/SoretDiffusionHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/modified/TrapMutationHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/modified/SoretDiffusionHandler.cpp
)
