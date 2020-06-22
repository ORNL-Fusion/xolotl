#ifndef W111MATERIALHANDLERFACTORY_H
#define W111MATERIALHANDLERFACTORY_H

#include <memory>

#include <xolotl/core/advection/W111AdvectionHandler.h>
#include <xolotl/core/flux/W111FitFluxHandler.h>
#include <xolotl/core/modified/W111TrapMutationHandler.h>
#include <xolotl/factory/material/MaterialFactory.h>

namespace xolotl
{
namespace factory
{
namespace material
{
/**
 * Subclass of MaterialFactory for a (111) oriented tungsten material.
 */
class W111MaterialFactory : public MaterialFactory
{
public:
	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W111MaterialFactory(const options::Options& options) :
		MaterialFactory(options)
	{
		theFluxHandler = std::make_shared<core::flux::W111FitFluxHandler>();
		theAdvectionHandler.push_back(
			std::make_shared<core::advection::W111AdvectionHandler>());
		theTrapMutationHandler =
			std::make_shared<core::modified::W111TrapMutationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	~W111MaterialFactory()
	{
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // W111MATERIALHANDLERFACTORY_H
