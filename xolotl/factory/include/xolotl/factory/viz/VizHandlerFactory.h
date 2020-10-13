#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace viz
{
class IVizHandler;
}

namespace factory
{
namespace viz
{
class VizHandlerFactory :
	public Factory<VizHandlerFactory, xolotl::viz::IVizHandler>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "VizHandlerFactory";
	}

	static std::string
	getName(const options::IOptions& options)
	{
		return options.getVizHandlerName();
	}
};
} // end namespace viz
} // end namespace factory
} // end namespace xolotl
