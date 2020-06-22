#ifndef W211MATERIALHANDLERFACTORY_H
#define W211MATERIALHANDLERFACTORY_H

#include <memory>

#include <xolotl/core/advection/W211AdvectionHandler.h>
#include <xolotl/core/flux/W211FitFluxHandler.h>
#include <xolotl/core/modified/W211TrapMutationHandler.h>
#include <xolotl/factory/material/MaterialFactory.h>

namespace xolotl
{
namespace factory
{
namespace material
{
/**
 * Subclass of MaterialFactory for a (211) oriented tungsten material.
 */
class W211MaterialFactory : public MaterialFactory
{
public:
	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W211MaterialFactory(const options::Options& options) :
		MaterialFactory(options)
	{
		theFluxHandler = std::make_shared<core::flux::W211FitFluxHandler>();
		theAdvectionHandler.push_back(
			std::make_shared<core::advection::W211AdvectionHandler>());
		theTrapMutationHandler =
			std::make_shared<core::modified::W211TrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W211MaterialFactory()
	{
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // W211MATERIALHANDLERFACTORY_H
