#pragma once

#include <xolotl/core/flux/VFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class VMaterialHandler : public MaterialHandler
{
public:
	VMaterialHandler(const options::IOptions& options) :
		MaterialHandler(
			options, MaterialSubHandlerGenerator<core::flux::VFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
