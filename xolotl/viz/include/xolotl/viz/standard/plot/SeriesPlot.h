#ifndef SERIESPLOT_H
#define SERIESPLOT_H

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
 * Plot the different data values as a function of one dimension. Available
 * PlottingStyle are POINTS or LINE. It can be associated to CvsXDataProvider.
 */
class SeriesPlot : public Plot
{
private:
	/**
	 * Container of data providers used for the plot.
	 */
	std::shared_ptr<std::vector<std::shared_ptr<dataprovider::IDataProvider>>>
		plotDataProviders;

public:
	/**
	 * The default constructor
	 */
	SeriesPlot() :
		Plot(),
		plotDataProviders(std::make_shared<
			std::vector<std::shared_ptr<dataprovider::IDataProvider>>>())
	{
	}

	/**
	 * The destructor
	 */
	~SeriesPlot()
	{
	}

	/**
	 * \see IPlot.h
	 */
	void
	render(const std::string& fileName = "fileName");

	/**
	 * \see IPlot.h
	 */
	void
	addDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider);

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider() const;

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider(int i) const;

	/**
	 * \see IPlot.h
	 */
	int
	getDataProviderNumber() const;
};

// end class SeriesPlot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
