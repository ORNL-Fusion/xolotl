#pragma once

#include <xolotl/core/advection/W211AdvectionHandler.h>
#include <xolotl/core/flux/W211FitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class W211MaterialHandler : public MaterialHandler
{
public:
	W211MaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::W211FitFluxHandler,
				core::advection::W211AdvectionHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
