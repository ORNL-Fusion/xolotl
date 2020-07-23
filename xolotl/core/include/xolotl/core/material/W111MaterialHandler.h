#pragma once

#include <xolotl/core/advection/W111AdvectionHandler.h>
#include <xolotl/core/flux/W111FitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/core/modified/W111TrapMutationHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class W111MaterialHandler : public MaterialHandler
{
public:
	W111MaterialHandler(const options::Options& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::W111FitFluxHandler,
				core::advection::W111AdvectionHandler,
				core::modified::W111TrapMutationHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
