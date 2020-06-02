#ifndef PLOTTINGSTYLE_H
#define PLOTTINGSTYLE_H

namespace xolotl {
namespace viz {

/**
 * Enumeration of the different plotting styles available.
 */
enum PlottingStyle {

	/**
	 * The data will be shown as points, whether it is a ScatterPlot, a SurfacePlot, or a VideoPlot.
	 */
	POINTS,

	/**
	 * Lines between the data points will be drawn for a ScatterPlot.
	 * A grid between the data points will be represented for a SurfacePlot or a VideoPlot.
	 */
	LINE,

	/**
	 * Only available for SurfacePlot and VideoPlot.
	 * The value at each data point will be represented by a color.
	 */
	COLORMAP,

	/**
	 * Only available for SurfacePlot and VideoPlot. The rendering will be a 3D surface as well as colors.
	 */
	GRID

};

//end enum PlottingStyle

} /* namespace viz */
} /* namespace xolotl */

#endif
