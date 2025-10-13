#pragma once

#include <xolotl/core/flux/RPVFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class RPVMaterialHandler : public MaterialHandler
{
public:
	RPVMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::RPVFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
