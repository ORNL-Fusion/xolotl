#ifndef SCATTERPLOT_H
#define SCATTERPLOT_H

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
 * Plot the data value as a function of one dimension. Available PlottingStyle
 * are POINTS or LINE. It can be associated to QvsXDataProvider,
 * QvsYDataProvider, QvsZDataProvider, or QvsTimeDataProvider.
 */
class ScatterPlot : public Plot
{
public:
	/**
	 * The default constructor
	 */
	ScatterPlot() : Plot()
	{
	}

	/**
	 * The destructor
	 */
	~ScatterPlot()
	{
	}

	/**
	 * \see IPlot.h
	 */
	void
	render(const std::string& fileName = "fileName");
};

// end class ScatterPlot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
