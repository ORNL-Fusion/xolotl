#include <xolotl/core/temperature/ConstantHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
namespace detail
{
auto constantTemperatureHandlerRegistration =
	xolotl::factory::temperature::TemperatureHandlerFactory::Registration<
		ConstantHandler>("constant");
}

ConstantHandler::ConstantHandler(double constTemperature) :
	temperature(constTemperature)
{
}

ConstantHandler::ConstantHandler(const options::IOptions& options) :
	ConstantHandler(options.getTempParam())
{
}

ConstantHandler::~ConstantHandler()
{
}
} // namespace temperature
} // namespace core
} // namespace xolotl
