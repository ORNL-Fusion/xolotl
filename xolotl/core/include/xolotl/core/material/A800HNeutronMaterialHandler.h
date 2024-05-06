#pragma once

#include <xolotl/core/flux/A800HNeutronFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class A800HNeutronMaterialHandler : public MaterialHandler
{
public:
	A800HNeutronMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::A800HNeutronFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
