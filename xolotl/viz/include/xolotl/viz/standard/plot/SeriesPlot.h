#ifndef SERIESPLOT_H
#define SERIESPLOT_H

// Includes
#include <vector>
#include <xolotl/viz/standard/plot/Plot.h>

namespace xolotl {
namespace viz {
namespace standard {
namespace plot {

/**
 * Plot the different data values as a function of one dimension. Available PlottingStyle are POINTS or LINE.
 * It can be associated to CvsXDataProvider.
 */
class SeriesPlot: public Plot {

private:

	/**
	 * Declare the constructor as private to force the use of a name
	 */
	SeriesPlot() :
			Plot() {
	}

private:

	/**
	 * Container of data providers used for the plot.
	 */
	std::shared_ptr<std::vector<std::shared_ptr<dataprovider::IDataProvider> > > plotDataProviders;

public:

	/**
	 * The default constructor
	 */
	SeriesPlot(const std::string& name);

	/**
	 * The destructor
	 */
	~SeriesPlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	void render(const std::string& fileName = "fileName");

	/**
	 * Method adding one data provider to the vector plotDataProviders
	 * @ param dataProvider The data provider to add.
	 */
	void addDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider);

	/**
	 * Method getting the i-th data provider
	 * @ param i The number of the data provider to be returned.
	 * @ return The ith data provider.
	 */
	std::shared_ptr<dataprovider::IDataProvider> getDataProvider(int i) const;

	/**
	 * Method getting the total number of data providers
	 * @ return The total number of data providers.
	 */
	int getDataProviderNumber() const;

};

//end class SeriesPlot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
