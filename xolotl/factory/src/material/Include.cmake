list(APPEND XOLOTL_FACTORY_HEADERS
    ${XOLOTL_FACTORY_HEADER_DIR}/material/AlloyMaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/FeMaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/FuelMaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/IMaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/MaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/PulsedMaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/TRIDYNMaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/W100MaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/W110MaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/W111MaterialFactory.h
    ${XOLOTL_FACTORY_HEADER_DIR}/material/W211MaterialFactory.h
)

list(APPEND XOLOTL_FACTORY_SOURCES
    ${XOLOTL_FACTORY_SOURCE_DIR}/material/IMaterialFactory.cpp
)
