#pragma once

#include <xolotl/core/flux/A800H9MeVFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class A800H9MeVMaterialHandler : public MaterialHandler
{
public:
	A800H9MeVMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::A800H9MeVFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
