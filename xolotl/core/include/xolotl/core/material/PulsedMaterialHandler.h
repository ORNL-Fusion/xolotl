#pragma once

#include <xolotl/core/flux/PulsedFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class PulsedMaterialHandler : public MaterialHandler
{
public:
	PulsedMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::PulsedFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
