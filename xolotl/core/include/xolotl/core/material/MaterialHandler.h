#pragma once

#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/flux/CustomFitFluxHandler.h>
#include <xolotl/core/material/IMaterialHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace material
{
/**
 * @brief Helper class interface for generating material sub-handlers for flux,
 * advection, and trap mutation
 *
 * @note We use this separate helper class (passed to the MaterialHandler
 * constructor) to avoid making MaterialHandler itself a class template.
 */
class IMaterialSubHandlerGenerator
{
public:
	/**
	 * @brief Generate the flux handler
	 *
	 * Expects flux handler constructor to take options
	 */
	virtual std::shared_ptr<core::flux::IFluxHandler>
	generateFluxHandler(const options::IOptions& options) const = 0;

	/**
	 * Generate the advection handler
	 */
	virtual std::shared_ptr<core::advection::IAdvectionHandler>
	generateAdvectionHandler() const = 0;

	/**
	 * Generate the trap mutation handler
	 */
	virtual std::shared_ptr<core::modified::ITrapMutationHandler>
	generateTrapMutationHandler() const = 0;
};

/**
 * @brief Helper class for generating the sub-handlers appropriate for a given
 * material
 *
 * MaterialHandler child classes must provide an instantiation of this template
 * in their constructors.
 *
 * @tparam TFluxHandler Must be a child of core::flux::IFluxHandler
 * @tparam TAdvectionHandler Must be a child of
 * core::advection::IAdvectionHandler
 * @tparam TTrapMutationHandler Must be a child of
 * core::modified::ITrapMutationHandler
 *
 * @note The flux handler type must be provided. However, if the options specify
 * a flux depth profile file, the CustomFitFluxHandler is used instead of the
 * specified flux handler type.
 */
template <typename TFluxHandler,
	typename TAdvectionHandler = core::advection::DummyAdvectionHandler,
	typename TTrapMutationHandler = core::modified::DummyTrapMutationHandler>
class MaterialSubHandlerGenerator : public IMaterialSubHandlerGenerator
{
	std::shared_ptr<core::flux::IFluxHandler>
	generateFluxHandler(const options::IOptions& options) const final
	{
		static_assert(
			std::is_base_of<core::flux::IFluxHandler, TFluxHandler>::value);

		if (options.getFluxDepthProfileFilePath().empty()) {
			return std::make_shared<TFluxHandler>(options);
		}
		else {
			return std::make_shared<flux::CustomFitFluxHandler>(options);
		}
	}

	std::shared_ptr<core::advection::IAdvectionHandler>
	generateAdvectionHandler() const final
	{
		static_assert(std::is_base_of<core::advection::IAdvectionHandler,
			TAdvectionHandler>::value);

		return std::make_shared<TAdvectionHandler>();
	}

	std::shared_ptr<core::modified::ITrapMutationHandler>
	generateTrapMutationHandler() const final
	{
		static_assert(std::is_base_of<core::modified::ITrapMutationHandler,
			TTrapMutationHandler>::value);

		return std::make_shared<TTrapMutationHandler>();
	}
};

/**
 * @brief MaterialHandler is a sort-of meta-handler class in that it enables a
 * set of handlers (flux, advection, diffusion, and trap mutation) to be defined
 * together for a given material
 */
class MaterialHandler : public IMaterialHandler
{
public:
	std::shared_ptr<core::flux::IFluxHandler>
	getFluxHandler() const final
	{
		return _fluxHandler;
	}

	std::vector<std::shared_ptr<core::advection::IAdvectionHandler>>
	getAdvectionHandler() const final
	{
		return _advectionHandlers;
	}

	std::shared_ptr<core::diffusion::IDiffusionHandler>
	getDiffusionHandler() const final
	{
		return _diffusionHandler;
	}

	std::shared_ptr<core::modified::ITrapMutationHandler>
	getTrapMutationHandler() const final
	{
		return _trapMutationHandler;
	}

protected:
	/**
	 * @brief Construct with options and sub-handler generator
	 *
	 * This constructor is intended to be used from a child class constructor.
	 * That is, a material handler implementation should provide the
	 * instantiation of MaterialSubHandlerGenerator to generate the appropriate
	 * sub-handler types.
	 */
	MaterialHandler(const options::IOptions& options,
		const IMaterialSubHandlerGenerator& subHandlerGenerator);

private:
	/**
	 * Generate diffusion handler appropriate for grid dimension if enabled
	 */
	static std::shared_ptr<core::diffusion::IDiffusionHandler>
	createDiffusionHandler(const options::IOptions& options);

	/**
	 * Finish setting up trap mutation handler based on relevant options
	 */
	void
	initializeTrapMutationHandler(const options::IOptions& options);

	/**
	 * Finish setting up advection handler(s) based on relevant options
	 */
	void
	initializeAdvectionHandlers(const options::IOptions& options);

protected:
	std::shared_ptr<core::diffusion::IDiffusionHandler> _diffusionHandler;
	std::vector<std::shared_ptr<core::advection::IAdvectionHandler>>
		_advectionHandlers;
	std::shared_ptr<core::flux::IFluxHandler> _fluxHandler;
	std::shared_ptr<core::modified::ITrapMutationHandler> _trapMutationHandler;
};
} // namespace material
} // namespace core
} // namespace xolotl
