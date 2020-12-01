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
	Plot() : plotStyle(PlottingStyle::LINE)
	{
	}

	/**
	 * The destructor.
	 */
	~Plot()
	{
	}

	/**
	 * \see IPlot.h
	 */
	void
	render(const std::string& fileName = "fileName")
	{
		return;
	}

	/**
	 * \see IPlot.h
	 */
	void
	write(const std::string& fileName)
	{
		return;
	}

	/**
	 * \see IPlot.h
	 */
	void
	setPlottingStyle(PlottingStyle style)
	{
		plotStyle = style;
		return;
	}

	/**
	 * \see IPlot.h
	 */
	PlottingStyle
	getPlottingStyle()
	{
		return plotStyle;
	}

	/**
	 * \see IPlot.h
	 */
	void
	setDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider)
	{
		plotDataProvider = dataProvider;
		return;
	}

	/**
	 * \see IPlot.h
	 */
	void
	addDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider)
	{
		plotDataProvider = dataProvider;
		return;
	}

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider() const
	{
		return plotDataProvider;
	}

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider(int i) const
	{
		return plotDataProvider;
	}

	/**
	 * \see IPlot.h
	 */
	int
	getDataProviderNumber() const
	{
		return 1;
	}

	/**
	 * \see IPlot.h
	 */
	void
	setLabelProvider(std::shared_ptr<LabelProvider> labelProvider)
	{
		plotLabelProvider = labelProvider;
		return;
	}

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<LabelProvider>
	getLabelProvider() const
	{
		return plotLabelProvider;
	}

	/**
	 * \see IPlot.h
	 */
	void
	showLegend(bool legendShow = true)
	{
		enableLegend = legendShow;
		return;
	}

	/**
	 * \see IPlot.h
	 */
	void
	setLogScale(bool logScale = true)
	{
		enableLogScale = logScale;
		return;
	}
};

// end class Plot

} /* namespace plot */
} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */

#endif
