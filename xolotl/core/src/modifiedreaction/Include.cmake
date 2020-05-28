list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/heterogeneousnucleation/DummyNucleationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/heterogeneousnucleation/HeterogeneousNucleationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/heterogeneousnucleation/IHeterogeneousNucleationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/DummyTrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/ITrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/Sigma3TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/W100TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/W110TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/W111TrapMutationHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/modifiedreaction/trapmutation/W211TrapMutationHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/modifiedreaction/heterogeneousnucleation/HeterogeneousNucleationHandler.cpp
    ${XOLOTL_CORE_SOURCE_DIR}/modifiedreaction/trapmutation/TrapMutationHandler.cpp
)
