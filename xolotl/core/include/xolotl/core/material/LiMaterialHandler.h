#pragma once

#include <xolotl/core/flux/LiFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class LiMaterialHandler : public MaterialHandler
{
public:
	LiMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::LiFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
