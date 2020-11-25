#ifndef DUMMYPLOT_H
#define DUMMYPLOT_H

#include <string>

#include <xolotl/viz/IPlot.h>

namespace xolotl
{
namespace viz
{
namespace dummy
{
/**
 * The DummyPlot class is instantiated by the DummyHandlerRegistry
 * class and realizes the IPlot interface.
 */
class DummyPlot : public IPlot
{
private:
	/**
	 * Declare the constructor as private to force the use of a name.
	 */
	DummyPlot()
	{
	}

	/**
	 * Data provider used for the plot.
	 */
	std::shared_ptr<dataprovider::IDataProvider> plotDataProvider;

public:
	/**
	 * The default constructor
	 */
	DummyPlot(const std::string& name);

	/**
	 * The destructor.
	 */
	~DummyPlot();

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 * \see IPlot.h
	 */
	void
	render(const std::string& fileName = "fileName");

	/**
	 * Method that will save the plotted plot in a file.
	 * \see IPlot.h
	 */
	void
	write(const std::string& fileName);

	/**
	 * Method allowing the user to set the PlottingStyle.
	 * \see IPlot.h
	 */
	void
	setPlottingStyle(PlottingStyle style);

	/**
	 * Method getting the PlottingStyle.
	 * \see IPlot.h
	 */
	PlottingStyle
	getPlottingStyle();

	/**
	 * Sets the data provider used for the plots.
	 * \see IPlot.h
	 */
	void
	setDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider);

	/**
	 * Method adding one data provider to the vector plotDataProviders
	 * \see IPlot.h
	 */
	void
	addDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider);

	/**
	 * Gets the data provider used.
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider() const;

	/**
	 * Method getting the i-th data provider for SeriesPlot
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider(int i) const;

	/**
	 * Method getting the total number of data providers
	 * \see IPlot.h
	 */
	int
	getDataProviderNumber() const;

	/**
	 * Sets the label provider used for the plots.
	 * \see IPlot.h
	 */
	void
	setLabelProvider(std::shared_ptr<LabelProvider> labelProvider);

	/**
	 * Gets the label provider used.
	 * \see IPlot.h
	 */
	std::shared_ptr<LabelProvider>
	getLabelProvider() const;

	/**
	 * Method that enables the rendering of the legend.
	 * \see IPlot.h
	 */
	void
	showLegend(bool legendShow = true);

	/**
	 * Method getting the legend.
	 * \see IPlot.h
	 */
	std::string
	getLegend() const;

	/**
	 * Method that enables the log scale.
	 * \see IPlot.h
	 */
	void
	setLogScale(bool logScale = true);
};
// end class DummyPlot

} // end namespace dummy
} // end namespace viz
} // end namespace xolotl

#endif
