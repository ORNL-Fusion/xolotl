#pragma once

#include <xolotl/core/flux/FeFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class FeMaterialHandler : public MaterialHandler
{
public:
	FeMaterialHandler(const options::Options& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::FeFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
