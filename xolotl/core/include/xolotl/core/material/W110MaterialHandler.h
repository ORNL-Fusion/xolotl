#pragma once

#include <xolotl/core/advection/W110AdvectionHandler.h>
#include <xolotl/core/flux/W110FitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class W110MaterialHandler : public MaterialHandler
{
public:
	W110MaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::W110FitFluxHandler,
				core::advection::W110AdvectionHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
