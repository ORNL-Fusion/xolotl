#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>

namespace xolotl
{
namespace factory
{
using core::temperature::ITemperatureHandler;
using temperature::TemperatureHandlerFactory;

template TemperatureHandlerFactory&
Factory<TemperatureHandlerFactory, ITemperatureHandler>::get();

template std::shared_ptr<ITemperatureHandler>
Factory<TemperatureHandlerFactory, ITemperatureHandler>::generate(
	const options::Options&);

template bool
Factory<TemperatureHandlerFactory, ITemperatureHandler>::registerGenerator(
	const std::string&, const Generator&);
} // end namespace factory
} // end namespace xolotl
