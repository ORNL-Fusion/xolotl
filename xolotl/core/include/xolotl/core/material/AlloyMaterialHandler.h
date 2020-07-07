#pragma once

#include <xolotl/core/flux/AlloyFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class AlloyMaterialHandler : public MaterialHandler
{
public:
	AlloyMaterialHandler(const options::Options& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::AlloyFitFluxHandler>{})
	{
	}
};

namespace detail
{
auto alloyMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		AlloyMaterialHandler>({"800H"});
}
} // namespace material
} // namespace core
} // namespace xolotl
