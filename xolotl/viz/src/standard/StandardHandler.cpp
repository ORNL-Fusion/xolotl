#include <stdexcept>

#include <xolotl/factory/viz/VizHandlerFactory.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/viz/LabelProvider.h>
#include <xolotl/viz/dataprovider/DataProvider.h>
#include <xolotl/viz/standard/StandardHandler.h>
#include <xolotl/viz/standard/plot/Plot.h>
#include <xolotl/viz/standard/plot/ScatterPlot.h>
#include <xolotl/viz/standard/plot/SeriesPlot.h>
#include <xolotl/viz/standard/plot/SurfacePlot.h>

namespace xolotl
{
namespace viz
{
namespace standard
{
#ifndef HAVE_VIZLIB_STD
namespace detail
{
auto stdHandlerRegistrations =
	xolotl::factory::viz::VizHandlerFactory::RegistrationCollection<
		StandardHandler>({"std"});
}
#endif

StandardHandler::StandardHandler(const options::IOptions& options)
{
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);
	if (procId == 0) {
		XOLOTL_LOG << "VizHandler: Using the standard handler";
	}
}

StandardHandler::~StandardHandler() = default;

std::shared_ptr<IPlot>
StandardHandler::getPlot(PlotType type)
{
	switch (type) {
	case PlotType::SCATTER:
		return std::make_shared<plot::ScatterPlot>();
	case PlotType::SERIES:
		return std::make_shared<plot::SeriesPlot>();
	case PlotType::SURFACE:
		return std::make_shared<plot::SurfacePlot>();
	default:
		return std::make_shared<plot::Plot>();
	}
}

} // namespace standard
} // namespace viz
} // namespace xolotl
