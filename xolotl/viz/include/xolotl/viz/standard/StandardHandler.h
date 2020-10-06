#pragma once

#include <string>

#include <xolotl/options/IOptions.h>
#include <xolotl/viz/IVizHandler.h>

namespace xolotl
{
namespace viz
{
namespace standard
{
/**
 * Handler for creating standard plots using EAVL and MESA libraries.
 * This is used only if the libraries are present and if the user uses
 * the standard registry.
 */
class StandardHandler : public IVizHandler
{
public:
	StandardHandler() = delete;

	StandardHandler(const options::IOptions& options);

	/**
	 * Clean up a StandardHandler.
	 */
	virtual ~StandardHandler();

	/**
	 * Obtain a Plot by name.
	 *
	 * @param name The name of the Plot.
	 * @param type The type of plot to return.
	 * @return A shared pointer to the newly-created Plot.
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(const std::string& name, PlotType type);
};
// end class StandardHandler

} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */
