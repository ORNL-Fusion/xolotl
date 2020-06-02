#ifndef VIDEOPLOT_H
#define VIDEOPLOT_H

// Includes
#include <vector>
#include <xolotl/viz/standard/plot/Plot.h>
#include <xolotl/core/Identifiable.h>

namespace xolotl {
namespace viz {
namespace standard {
namespace plot {

/**
 * Plot the data value as a function of two different spatial dimensions for each time step
 * and change the time step with time to have a video-like rendering, each frame being a SurfacePlot.
 * The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE.
 * It can be associated to QvsXYTimeDataProvider, QvsXZTimeDataProvider, or QvsYZTimeDataProvider.
 */
class VideoPlot: public Plot {

private:

	/**
	 * Declare the constructor as private to force the use of a name
	 */
	VideoPlot() :
			Plot() {
	}

private:

	/**
	 * Number of frames shown per second.
	 */
	double frameRate;

public:

	/**
	 * The default constructor
	 */
	VideoPlot(const std::string& name);

	/**
	 * The destructor
	 */
	~VideoPlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	void render(const std::string& fileName = "fileName");

	/**
	 * Sets the frame rate for VideoPlot.
	 */
	void setFrameRate(double frameRate);

	/**
	 * Get the frame rate from VideoPlot.
	 */
	double getFrameRate() const;

};

//end class VideoPlot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
