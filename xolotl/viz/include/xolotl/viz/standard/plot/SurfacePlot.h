#ifndef SURFACEPLOT_H
#define SURFACEPLOT_H

// Includes
#include <vector>

#include <xolotl/viz/standard/plot/Plot.h>

namespace xolotl
{
namespace viz
{
namespace standard
{
namespace plot
{
/**
 * Plot the data value as a function of two different dimensions.
 * The available PlottingStyle are POINTS, LINE, COLORMAP, and SURFACE.
 * It can be associated to CvsXYDataProvider.
 */
class SurfacePlot : public Plot
{
public:
	/**
	 * The default constructor
	 */
	SurfacePlot() : Plot()
	{
	}

	/**
	 * The destructor
	 */
	~SurfacePlot()
	{
	}

	/**
	 * \see IPlot.h
	 */
	void
	render(const std::string& fileName = "fileName");
};

// end class SurfacePlot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
