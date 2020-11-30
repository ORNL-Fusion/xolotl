#pragma once

#include <xolotl/core/flux/FuelFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class FuelMaterialHandler : public MaterialHandler
{
public:
	FuelMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::FuelFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
