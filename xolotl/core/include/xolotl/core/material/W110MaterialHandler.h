#pragma once

#include <xolotl/core/advection/W110AdvectionHandler.h>
#include <xolotl/core/flux/W110FitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/core/modified/W110TrapMutationHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class W110MaterialHandler : public MaterialHandler
{
public:
	W110MaterialHandler(const options::Options& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::W110FitFluxHandler,
				core::advection::W110AdvectionHandler,
				core::modified::W110TrapMutationHandler>{})
	{
	}
};

namespace detail
{
auto w110MaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		W110MaterialHandler>({"W110"});
}
} // namespace material
} // namespace core
} // namespace xolotl
