#ifndef FUELMATERIALHANDLERFACTORY_H
#define FUELMATERIALHANDLERFACTORY_H

#include <memory>

#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/flux/FuelFitFluxHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>
#include <xolotl/factory/material/MaterialFactory.h>

namespace xolotl
{
namespace factory
{
namespace material
{
/**
 * Subclass of MaterialFactory for a nuclear fuel material.
 */
class FuelMaterialFactory : public MaterialFactory
{
public:
	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	FuelMaterialFactory(const options::Options& options) :
		MaterialFactory(options)
	{
		theFluxHandler = std::make_shared<core::flux::FuelFitFluxHandler>();
		theAdvectionHandler.push_back(
			std::make_shared<core::advection::DummyAdvectionHandler>());
		theTrapMutationHandler =
			std::make_shared<core::modified::DummyTrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~FuelMaterialFactory()
	{
	}

	/**
	 * Initialize the material conditions with the different given options.
	 * \see IMaterialFactory.h
	 */
	void
	initializeMaterial(const options::Options& options)
	{
		// Call the mother method
		MaterialFactory::initializeMaterial(options);

		// Change the flux amplitude because we have to take into account
		// that there are one xenon created every 4 fissions.
		theFluxHandler->setFluxAmplitude(
			options.getFluxAmplitude() * options.getFissionYield());

		return;
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // FUELMATERIALHANDLERFACTORY_H
