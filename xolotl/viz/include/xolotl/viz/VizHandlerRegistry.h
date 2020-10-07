#pragma once

#include <memory>

#include <xolotl/viz/IVizHandler.h>

namespace xolotl
{
namespace viz
{
class VizHandlerRegistry
{
public:
	static void
	set(const std::shared_ptr<IVizHandler>& handler);

	static IVizHandler*
	get()
	{
		return vizHandler.get();
	}

private:
	static std::shared_ptr<IVizHandler> vizHandler;
};
} // namespace viz
} // namespace xolotl
