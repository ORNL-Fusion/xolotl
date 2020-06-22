// Includes
#include <xolotl/viz/dummy/DummyPlot.h>

namespace xolotl
{
namespace viz
{
namespace dummy
{
DummyPlot::DummyPlot(const std::string&) : util::Identifiable("unused")
{
}

DummyPlot::~DummyPlot()
{
}

void
DummyPlot::render(const std::string&)
{
	return;
}

void
DummyPlot::write(const std::string&)
{
	return;
}

void DummyPlot::setPlottingStyle(PlottingStyle)
{
	return;
}

PlottingStyle
DummyPlot::getPlottingStyle()
{
	return PlottingStyle();
}

void
DummyPlot::setDataProvider(
	std::shared_ptr<dataprovider::IDataProvider> dataProvider)
{
	plotDataProvider = dataProvider;
	return;
}

std::shared_ptr<dataprovider::IDataProvider>
DummyPlot::getDataProvider() const
{
	return plotDataProvider;
}

void
DummyPlot::setLabelProvider(std::shared_ptr<LabelProvider> labelProvider)
{
	plotLabelProvider = labelProvider;
	return;
}

std::shared_ptr<LabelProvider>
DummyPlot::getLabelProvider() const
{
	return plotLabelProvider;
}

void
DummyPlot::showLegend(bool)
{
	return;
}

std::string
DummyPlot::getLegend() const
{
	return " ";
}

void
DummyPlot::setLogScale(bool)
{
	return;
}

void
DummyPlot::addDataProvider(
	std::shared_ptr<dataprovider::IDataProvider> dataProvider)
{
	plotDataProvider = dataProvider;
}

std::shared_ptr<dataprovider::IDataProvider>
DummyPlot::getDataProvider(int) const
{
	return plotDataProvider;
}

int
DummyPlot::getDataProviderNumber() const
{
	return 0;
}

} // namespace dummy
} // namespace viz
} // namespace xolotl
