#pragma once

#include <xolotl/core/advection/W100AdvectionHandler.h>
#include <xolotl/core/flux/W100FitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/core/modified/SoretDiffusionHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class W100MaterialHandler : public MaterialHandler
{
public:
	W100MaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::W100FitFluxHandler,
				core::advection::W100AdvectionHandler,
				core::modified::SoretDiffusionHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
