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
class T91MaterialHandler : public MaterialHandler
{
public:
	T91MaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::FeFitFluxHandler>{})
	{
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
