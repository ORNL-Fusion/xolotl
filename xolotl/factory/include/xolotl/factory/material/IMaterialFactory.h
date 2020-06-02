#ifndef IMATERIALHANDLERFACTORY_H
#define IMATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/options/Options.h>
#include <xolotl/core/flux/IFluxHandler.h>
#include <xolotl/core/advection/IAdvectionHandler.h>
#include <xolotl/core/diffusion/IDiffusionHandler.h>
#include <xolotl/core/modified/ITrapMutationHandler.h>
#include <xolotl/core/modified/IHeterogeneousNucleationHandler.h>

namespace xolotl {
namespace factory {
namespace material {

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
	virtual void initializeMaterial(const options::Options &opts) = 0;

	/**
	 * Return the flux handler.
	 *
	 * @return The flux handler.
	 */
	virtual std::shared_ptr<core::flux::IFluxHandler> getFluxHandler() const = 0;

	/**
	 * Return the advection handlers.
	 *
	 * @return The advection handlers.
	 */
	virtual std::vector<std::shared_ptr<core::advection::IAdvectionHandler> > getAdvectionHandler() const = 0;

	/**
	 * Return the diffusion handler.
	 *
	 * @return The diffusion handler.
	 */
	virtual std::shared_ptr<core::diffusion::IDiffusionHandler> getDiffusionHandler() const = 0;

	/**
	 * Return the modified trap-mutation handler.
	 *
	 * @return The trap mutation handler.
	 */
	virtual std::shared_ptr<core::modified::ITrapMutationHandler> getTrapMutationHandler() const = 0;

	/**
	 * Return the heterogeneous nucleation handler.
	 *
	 * @return The nucleation handler.
	 */
	virtual std::shared_ptr<core::modified::IHeterogeneousNucleationHandler> getNucleationHandler() const = 0;

	/**
	 * Function that create the wanted material factory depending on the given type.
	 *
	 * @param options The options
	 * @return The material factory.
	 */
	static std::shared_ptr<IMaterialFactory> createMaterialFactory(
			const options::Options &options);

};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // IMATERIALHANDLERFACTORY_H
