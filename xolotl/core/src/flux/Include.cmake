list(APPEND XOLOTL_CORE_HEADERS
    ${XOLOTL_CORE_HEADER_DIR}/flux/AlloyFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/AlloySRIMData.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/AlphaZrFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/CustomFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/FeFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/FeCrFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/FluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/FuelFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/IFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/PSIFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/PulsedFitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/W100FitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/W110FitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/W111FitFluxHandler.h
    ${XOLOTL_CORE_HEADER_DIR}/flux/W211FitFluxHandler.h
)

list(APPEND XOLOTL_CORE_SOURCES
    ${XOLOTL_CORE_SOURCE_DIR}/flux/FluxHandler.cpp
)
