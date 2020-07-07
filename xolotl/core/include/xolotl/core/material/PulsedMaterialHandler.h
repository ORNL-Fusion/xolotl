#pragma once

#include <xolotl/core/flux/PulsedFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class PulsedMaterialHandler : public MaterialHandler
{
public:
	PulsedMaterialHandler(const options::Options& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::PulsedFitFluxHandler>{})
	{
		_fluxHandler->setPulseTime(options.getPulseTime());
		_fluxHandler->setProportion(options.getPulseProportion());
	}
};

namespace detail
{
auto pulsedMaterialHandlerRegistrations =
	xolotl::factory::material::MaterialHandlerFactory::RegistrationCollection<
		PulsedMaterialHandler>({"Pulsed"});
}
} // namespace material
} // namespace core
} // namespace xolotl
