#include <xolotl/core/temperature/GradientHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
namespace detail
{
auto gradientTemperatureHandlerRegistration =
	xolotl::factory::temperature::TemperatureHandlerFactory::Registration<
		GradientHandler>("gradient");
}

GradientHandler::GradientHandler(double surfTemp, double bulkTemp) :
	surfaceTemperature(surfTemp),
	bulkTemperature(bulkTemp)
{
}

GradientHandler::GradientHandler(const options::IOptions& options) :
	GradientHandler(options.getTempParam(0), options.getTempParam(1))
{
}

GradientHandler::~GradientHandler()
{
}

double
GradientHandler::getTemperature(
	const plsm::SpaceVector<double, 3>& fraction, double) const
{
	return surfaceTemperature +
		(bulkTemperature - surfaceTemperature) * fraction[0];
}
} // namespace temperature
} // namespace core
} // namespace xolotl
