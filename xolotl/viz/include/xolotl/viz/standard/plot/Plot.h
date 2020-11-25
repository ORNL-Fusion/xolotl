#ifndef PLOT_H
#define PLOT_H

// Includes
#include <xolotl/viz/IPlot.h>

namespace xolotl
{
namespace viz
{
namespace standard
{
namespace plot
{
/**
 * Plot is the class that realizes the interface IPlot.
 * It is a general class that provides general methods, but to actual plot
 * anything, the user needs to use one of its subclasses.
 */
class Plot : public IPlot
{
protected:
	/**
	 * Declare the constructor as private to force the use of a name
	 */
	Plot()
	{
	}

	/**
	 * Choice of PlottingStyle.
	 */
	PlottingStyle plotStyle;

	/**
	 * If it is equal to True, the legend will be displayed.
	 */
	bool enableLegend = false;

	/**
	 * If it is equal to True, a log scale will be used (for 1D plot for now).
	 */
	bool enableLogScale = false;

	/**
	 * Data provider used for the plot.
	 */
	std::shared_ptr<dataprovider::IDataProvider> plotDataProvider;

public:
	/**
	 * The default constructor
	 */
	Plot(const std::string& name);

	/**
	 * The destructor.
	 */
	~Plot();

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

// end class Plot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
