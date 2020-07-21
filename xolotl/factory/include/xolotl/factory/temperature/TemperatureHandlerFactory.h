#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
class ITemperatureHandler;
}
} // namespace core

namespace factory
{
namespace temperature
{
class TemperatureHandlerFactory :
	public Factory<TemperatureHandlerFactory,
		core::temperature::ITemperatureHandler>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "TemperatureHandlerFactory";
	}

	static std::string
	getName(const options::Options&)
	{
		return "";
	}
};
} // end namespace temperature
} // end namespace factory
} // end namespace xolotl
