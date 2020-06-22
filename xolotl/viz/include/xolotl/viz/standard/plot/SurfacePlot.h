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
 * It can be associated to QvsXYDataProvider, QvsXZDataProvider, or
 * QvsYZDataProvider.
 */
class SurfacePlot : public Plot
{
private:
	/**
	 * Declare the constructor as private to force the use of a name
	 */
	SurfacePlot() : Plot()
	{
	}

public:
	/**
	 * The default constructor
	 */
	SurfacePlot(const std::string& name);

	/**
	 * The destructor
	 */
	~SurfacePlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
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
