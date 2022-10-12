#pragma once

#include <xolotl/core/flux/AlphaZrFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class AlphaZrMaterialHandler : public MaterialHandler
{
public:
	AlphaZrMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::AlphaZrFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
