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
public:
	/**
	 * The default constructor
	 */
	DummyPlot()
	{
	}

	/**
	 * The destructor.
	 */
	~DummyPlot()
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
		return;
	}

	/**
	 * \see IPlot.h
	 */
	PlottingStyle
	getPlottingStyle()
	{
		return PlottingStyle();
	}

	/**
	 * \see IPlot.h
	 */
	void
	setDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider)
	{
		return;
	}

	/**
	 * \see IPlot.h
	 */
	void
	addDataProvider(std::shared_ptr<dataprovider::IDataProvider> dataProvider)
	{
		return;
	}

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider() const
	{
		return std::shared_ptr<dataprovider::IDataProvider>();
	}

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<dataprovider::IDataProvider>
	getDataProvider(int i) const
	{
		return std::shared_ptr<dataprovider::IDataProvider>();
	}

	/**
	 * \see IPlot.h
	 */
	int
	getDataProviderNumber() const
	{
		return 0;
	}

	/**
	 * \see IPlot.h
	 */
	void
	setLabelProvider(std::shared_ptr<LabelProvider> labelProvider)
	{
		return;
	}

	/**
	 * \see IPlot.h
	 */
	std::shared_ptr<LabelProvider>
	getLabelProvider() const
	{
		return std::shared_ptr<LabelProvider>();
	}

	/**
	 * \see IPlot.h
	 */
	void
	showLegend(bool legendShow = true)
	{
		return;
	}

	/**
	 * \see IPlot.h
	 */
	void
	setLogScale(bool logScale = true)
	{
		return;
	}
};
// end class DummyPlot

} // end namespace dummy
} // end namespace viz
} // end namespace xolotl

#endif
