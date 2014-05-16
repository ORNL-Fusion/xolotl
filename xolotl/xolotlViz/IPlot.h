#ifndef IPLOT_H
#define IPLOT_H

// Includes
#include "PlottingStyle.h"
#include <LabelProvider.h>
#include <IDataProvider.h>
#include <string>
#include <memory>
#include <IIdentifiable.h>

namespace xolotlViz {

/**
 * IPlot describe the structure needed to be able to plot data provided by DataProvider.
 * The user interacts with it through different method where he/she could set the data to plot,
 * title, legend, plotting style, etc.
 */
class IPlot : public virtual xolotlCore::IIdentifiable {

public:

	/**
	 * LabelProvider used for the Plot.
	 */
	std::shared_ptr<LabelProvider> plotLabelProvider;

	/**
	 * The destructor
	 */
	virtual ~IPlot() { }

	/**
	 * Method managing everything that is related to the rendering of a plot.
	 */
	virtual void render(std::string fileName = "fileName") = 0;

	/**
	 * Method that will save the plotted plot in a file.
	 * @param fileName The name of the file where the plot
	 * will be saved.
	 */
	virtual void write(std::string fileName) = 0;

	/**
	 * Method allowing the user to set the PlottingStyle.
	 * @param style The PlottingStyle that will be used.
	 */
	virtual void setPlottingStyle(PlottingStyle style) = 0;

	/**
	 * Method getting the PlottingStyle.
	 * @return The PlottingStyle attribute.
	 */
	virtual PlottingStyle getPlottingStyle() = 0;

	/**
	 * Sets the data provider used for the plots.
	 * @param dataProvider The DataProvider to set.
	 */
	virtual void setDataProvider(std::shared_ptr<IDataProvider> dataProvider) = 0;

	/**
	 * Method adding one data provider to the vector plotDataProviders
	 * @ param dataProvider The data provider to add.
	 */
	virtual void addDataProvider(std::shared_ptr<IDataProvider> dataProvider) = 0;

	/**
	 * Gets the data provider used.
	 * @return the DataProvider attribute.
	 */
	virtual std::shared_ptr<IDataProvider> getDataProvider() const = 0;

	/**
	 * Method getting the i-th data provider for SeriesPlot
	 * @ param i The number of the data provider to be returned.
	 * @ return The ith data provider.
	 */
	virtual std::shared_ptr<IDataProvider> getDataProvider(int i) const = 0;

	/**
	 * Method getting the total number of data providers
	 * @ return The total number of data providers.
	 */
	virtual int getDataProviderNumber() const = 0;

	/**
	 * Sets the label provider used for the plots.
	 * @param dataProvider The LabelProvider to set.
	 */
	virtual void setLabelProvider(std::shared_ptr<LabelProvider> labelProvider) = 0;

	/**
	 * Gets the label provider used.
	 * @return the LabelProvider attribute.
	 */
	virtual std::shared_ptr<LabelProvider> getLabelProvider() const = 0;

	/**
	 * Method that enables the rendering of the legend.
	 * @param legendShow A boolean set to true by default
	 */
	virtual void showLegend(bool legendShow = true) = 0;

	/**
	 * Method defining the legend with the help of the data provider and the label provider.
	 * @return The legend.
	 */
	virtual std::string getLegend() const = 0;

	/**
	 * Method that enables the log scale.
	 * @param logScale A boolean set to true by default
	 */
	virtual void setLogScale(bool logScale = true) = 0;

};

//end class IPlot

} /* namespace xolotlViz */

#endif
