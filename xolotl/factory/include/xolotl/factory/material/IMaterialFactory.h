#ifndef IMATERIALHANDLERFACTORY_H
#define IMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/options/Options.h>
#include <xolotl/core/flux/IFluxHandler.h>
#include <xolotl/core/advection/IAdvectionHandler.h>
#include <xolotl/core/diffusion/IDiffusionHandler.h>
#include <xolotl/core/modifiedreaction/trapmutation/ITrapMutationHandler.h>
#include <xolotl/core/modifiedreaction/heterogeneousnucleation/IHeterogeneousNucleationHandler.h>

namespace xolotlFactory {

/**
 * Realizations of this interface are responsible for handling the flux and the advection.
 * they are both dependent on the type of material under study.
 */
class IMaterialFactory {
public:

	/**
	 * The destructor
	 */
	~IMaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	virtual void initializeMaterial(const xolotlCore::Options &options) = 0;

	/**
	 * Return the flux handler.
	 *
	 * @return The flux handler.
	 */
	virtual std::shared_ptr<xolotlCore::IFluxHandler> getFluxHandler() const = 0;

	/**
	 * Return the advection handlers.
	 *
	 * @return The advection handlers.
	 */
	virtual std::vector<std::shared_ptr<xolotlCore::IAdvectionHandler> > getAdvectionHandler() const = 0;

	/**
	 * Return the diffusion handler.
	 *
	 * @return The diffusion handler.
	 */
	virtual std::shared_ptr<xolotlCore::IDiffusionHandler> getDiffusionHandler() const = 0;

	/**
	 * Return the modified trap-mutation handler.
	 *
	 * @return The trap mutation handler.
	 */
	virtual std::shared_ptr<xolotlCore::ITrapMutationHandler> getTrapMutationHandler() const = 0;

	/**
	 * Return the heterogeneous nucleation handler.
	 *
	 * @return The nucleation handler.
	 */
	virtual std::shared_ptr<xolotlCore::IHeterogeneousNucleationHandler> getNucleationHandler() const = 0;

	/**
	 * Function that create the wanted material factory depending on the given type.
	 *
	 * @param options The options
	 * @return The material factory.
	 */
	static std::shared_ptr<IMaterialFactory> createMaterialFactory(
			const xolotlCore::Options &options);

};

} // end namespace xolotlFactory

#endif // IMATERIALHANDLERFACTORY_H
