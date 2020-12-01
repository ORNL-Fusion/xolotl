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
 * Handler for creating standard plots using VTK-m library.
 * This is used only if the library is present and if the user uses
 * the standard registry.
 */
class StandardHandler : public IVizHandler
{
public:
	StandardHandler() = delete;

	StandardHandler(const options::IOptions& options);

	virtual ~StandardHandler();

	/**
	 * \see IVizHandler.h
	 */
	virtual std::shared_ptr<IPlot>
	getPlot(PlotType type);
};
// end class StandardHandler

} /* namespace standard */
} /* namespace viz */
} /* namespace xolotl */
