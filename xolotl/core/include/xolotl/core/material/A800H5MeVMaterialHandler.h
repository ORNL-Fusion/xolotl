#pragma once

#include <xolotl/core/flux/A800H5MeVFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class A800H5MeVMaterialHandler : public MaterialHandler
{
public:
	A800H5MeVMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::A800H5MeVFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
