#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace viz
{
class IVizHandlerRegistry;
}

namespace factory
{
namespace viz
{
class VizHandlerRegistryFactory :
	public Factory<VizHandlerRegistryFactory, xolotl::viz::IVizHandlerRegistry>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "VizHandlerRegistryFactory";
	}

	static std::string
	getName(const options::Options& options)
	{
		return "";
	}
};
} // end namespace viz
} // end namespace factory
} // end namespace xolotl
