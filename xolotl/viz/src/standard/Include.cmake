list(APPEND XOLOTL_VIZ_HEADERS
    ${XOLOTL_VIZ_HEADER_DIR}/standard/StandardHandler.h
    ${XOLOTL_VIZ_HEADER_DIR}/standard/plot/Plot.h
    ${XOLOTL_VIZ_HEADER_DIR}/standard/plot/ScatterPlot.h
    ${XOLOTL_VIZ_HEADER_DIR}/standard/plot/SeriesPlot.h
    ${XOLOTL_VIZ_HEADER_DIR}/standard/plot/SurfacePlot.h
)

list(APPEND XOLOTL_VIZ_SOURCES
    ${XOLOTL_VIZ_SOURCE_DIR}/standard/StandardHandler.cpp
    ${XOLOTL_VIZ_SOURCE_DIR}/standard/plot/ScatterPlot.cpp
    ${XOLOTL_VIZ_SOURCE_DIR}/standard/plot/SeriesPlot.cpp
    ${XOLOTL_VIZ_SOURCE_DIR}/standard/plot/SurfacePlot.cpp
)
