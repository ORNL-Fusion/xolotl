#pragma once

#include <xolotl/core/flux/FeCrFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class FeCrMaterialHandler : public MaterialHandler
{
public:
	FeCrMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::FeCrFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
