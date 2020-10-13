#pragma once

#include <xolotl/core/flux/FuelFitFluxHandler.h>
#include <xolotl/core/material/MaterialHandler.h>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class FuelMaterialHandler : public MaterialHandler
{
public:
	FuelMaterialHandler(const options::IOptions& options) :
		MaterialHandler(options,
			MaterialSubHandlerGenerator<core::flux::FuelFitFluxHandler>{})
	{
		// Change the flux amplitude because we have to take into account
		// that there are one xenon created every 4 fissions.
		_fluxHandler->setFluxAmplitude(
			options.getFluxAmplitude() * options.getFissionYield());
	}
};
} // namespace material
} // namespace core
} // namespace xolotl
