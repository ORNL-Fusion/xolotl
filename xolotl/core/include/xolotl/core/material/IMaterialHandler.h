#pragma once

#include <memory>
#include <vector>

#include <xolotl/core/advection/IAdvectionHandler.h>
#include <xolotl/core/diffusion/IDiffusionHandler.h>
#include <xolotl/core/flux/IFluxHandler.h>
#include <xolotl/core/modified/ITrapMutationHandler.h>

namespace xolotl
{
namespace core
{
namespace material
{
class IMaterialHandler
{
public:
	virtual ~IMaterialHandler()
	{
	}

	/**
	 * Return the flux handler.
	 *
	 * @return The flux handler.
	 */
	virtual std::shared_ptr<core::flux::IFluxHandler>
	getFluxHandler() const = 0;

	/**
	 * Return the advection handlers.
	 *
	 * @return The advection handlers.
	 *
	 * @todo TODO: Should this have an index parameter and just return the
	 * specified handler?
	 */
	virtual std::vector<std::shared_ptr<core::advection::IAdvectionHandler>>
	getAdvectionHandler() const = 0;

	/**
	 * Return the diffusion handler.
	 *
	 * @return The diffusion handler.
	 */
	virtual std::shared_ptr<core::diffusion::IDiffusionHandler>
	getDiffusionHandler() const = 0;

	/**
	 * Return the modified trap-mutation handler.
	 *
	 * @return The trap mutation handler.
	 */
	virtual std::shared_ptr<core::modified::ITrapMutationHandler>
	getTrapMutationHandler() const = 0;
};
} // namespace material
} // namespace core
} // namespace xolotl
