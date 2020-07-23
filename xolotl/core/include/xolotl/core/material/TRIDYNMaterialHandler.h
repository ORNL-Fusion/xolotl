#pragma once

#include <xolotl/core/advection/W100AdvectionHandler.h>
#include <xolotl/core/flux/TRIDYNFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/core/modified/W100TrapMutationHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class TRIDYNMaterialHandler : public MaterialHandler
{
public:
	TRIDYNMaterialHandler(const options::Options& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::TRIDYNFitFluxHandler,
				core::advection::W100AdvectionHandler,
				core::modified::W100TrapMutationHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
