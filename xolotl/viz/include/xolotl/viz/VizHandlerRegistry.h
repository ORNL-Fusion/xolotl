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
	/**
	 * This operation sets the visualization handler.
	 *
	 * @param handler
	 */
	static void
	set(const std::shared_ptr<IVizHandler>& handler);

	/**
	 * This operation gets the visualization handler.
	 *
	 * @return handler
	 */
	static const std::shared_ptr<IVizHandler>&
	get();

private:
	static std::shared_ptr<IVizHandler> vizHandler;
};
} // namespace viz
} // namespace xolotl
