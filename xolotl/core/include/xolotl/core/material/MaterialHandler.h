#pragma once

#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/material/IMaterialHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>
#include <xolotl/options/Options.h>

namespace xolotl
{
namespace core
{
namespace material
{
class IMaterialSubHandlerGenerator
{
public:
	virtual std::shared_ptr<core::flux::IFluxHandler>
	generateFluxHandler() const = 0;

	virtual std::shared_ptr<core::advection::IAdvectionHandler>
	generateAdvectionHandler() const = 0;

	virtual std::shared_ptr<core::modified::ITrapMutationHandler>
	generateTrapMutationHandler() const = 0;
};

template <typename TFluxHandler,
	typename TAdvectionHandler = core::advection::DummyAdvectionHandler,
	typename TTrapMutationHandler = core::modified::DummyTrapMutationHandler>
class MaterialSubHandlerGenerator : public IMaterialSubHandlerGenerator
{
	std::shared_ptr<core::flux::IFluxHandler>
	generateFluxHandler() const final
	{
		static_assert(
			std::is_base_of<core::flux::IFluxHandler, TFluxHandler>::value, "");

		return std::make_shared<TFluxHandler>();
	}

	std::shared_ptr<core::advection::IAdvectionHandler>
	generateAdvectionHandler() const final
	{
		static_assert(std::is_base_of<core::advection::IAdvectionHandler,
						  TAdvectionHandler>::value,
			"");

		return std::make_shared<TAdvectionHandler>();
	}

	std::shared_ptr<core::modified::ITrapMutationHandler>
	generateTrapMutationHandler() const final
	{
		static_assert(std::is_base_of<core::modified::ITrapMutationHandler,
						  TTrapMutationHandler>::value,
			"");

		return std::make_shared<TTrapMutationHandler>();
	}
};

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
	MaterialHandler(const options::Options& options,
		const IMaterialSubHandlerGenerator& subHandlerGenerator);

private:
	static std::shared_ptr<core::diffusion::IDiffusionHandler>
	createDiffusionHandler(const options::Options& options);

	void
	initializeFluxHandler(const options::Options& options);

	void
	initializeTrapMutationHandler(const options::Options& options);

	void
	initializeAdvectionHandlers(const options::Options& options);

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
