#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include <xolotl/factory/material/IMaterialFactory.h>
#include <xolotl/core/diffusion/DummyDiffusionHandler.h>
#include <xolotl/core/advection/DummyAdvectionHandler.h>
#include <xolotl/core/modified/DummyTrapMutationHandler.h>
#include <xolotl/core/modified/DummyNucleationHandler.h>
#include <xolotl/core/advection/XGBAdvectionHandler.h>
#include <xolotl/core/advection/YGBAdvectionHandler.h>
#include <xolotl/core/advection/ZGBAdvectionHandler.h>
#include <xolotl/core/diffusion/Diffusion1DHandler.h>
#include <xolotl/core/diffusion/Diffusion2DHandler.h>
#include <xolotl/core/diffusion/Diffusion3DHandler.h>
#include <xolotl/util/TokenizedLineReader.h>

namespace xolotl {
namespace factory {
namespace material {

/**
 * Realizes the IMaterialFactory interface. Handles the flux, the advection, modified
 * trap mutation, and for a specific material.
 */
class MaterialFactory: public IMaterialFactory {
protected:

	//! The flux handler
	std::shared_ptr<core::flux::IFluxHandler> theFluxHandler;

	//! The advection handler
	std::vector<std::shared_ptr<core::advection::IAdvectionHandler> > theAdvectionHandler;

	//! The diffusion handler
	std::shared_ptr<core::diffusion::IDiffusionHandler> theDiffusionHandler;

	//! The modified trap-mutation handler
	std::shared_ptr<core::modified::ITrapMutationHandler> theTrapMutationHandler;

	//! The heterogeneous nucleation handler
	std::shared_ptr<core::modified::IHeterogeneousNucleationHandler> theNucleationHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	MaterialFactory(const options::Options &opts) {
		// Get the dimension and migration energy threshold
		int dim = opts.getDimensionNumber();
		double threshold = opts.getMigrationThreshold();
		// Switch on the dimension for the diffusion handler
		switch (dim) {
		case 0:
			theDiffusionHandler = std::make_shared<
					core::diffusion::DummyDiffusionHandler>(threshold);
			break;
		case 1:
			theDiffusionHandler = std::make_shared<
					core::diffusion::Diffusion1DHandler>(threshold);
			break;
		case 2:
			theDiffusionHandler = std::make_shared<
					core::diffusion::Diffusion2DHandler>(threshold);
			break;
		case 3:
			theDiffusionHandler = std::make_shared<
					core::diffusion::Diffusion3DHandler>(threshold);
			break;
		default:
			// The asked dimension is not good (e.g. -1, 4)
			throw std::string(
					"\nxolotlFactory: Bad dimension for the material factory.");
		}
	}

	/**
	 * The destructor
	 */
	~MaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	virtual void initializeMaterial(const options::Options &opts) {
		// Wrong if both he flux and time profile options are used
		if (opts.useFluxAmplitude() && opts.useFluxTimeProfile()) {
			// A constant flux value AND a time profile cannot both be given.
			throw std::string(
					"\nA constant flux value AND a time profile cannot both be given.");
		} else if (opts.useFluxAmplitude()) {
			// Set the constant value of the flux
			theFluxHandler->setFluxAmplitude(opts.getFluxAmplitude());
		} else if (opts.useFluxTimeProfile()) {
			// Initialize the time profile
			theFluxHandler->initializeTimeProfile(opts.getFluxProfileName());
		}

		// Get the process map
		auto map = opts.getProcesses();
		// Set dummy handlers when needed
		if (!map["diff"])
			theDiffusionHandler = std::make_shared<
					core::diffusion::DummyDiffusionHandler>(
					opts.getMigrationThreshold());
		if (!map["advec"]) {
			// Clear the advection handler
			theAdvectionHandler.clear();
			// To replace it by a dummy one
			theAdvectionHandler.push_back(
					std::make_shared<core::advection::DummyAdvectionHandler>());
		}
		if (!map["modifiedTM"])
			theTrapMutationHandler = std::make_shared<
					core::modified::DummyTrapMutationHandler>();
		if (!map["attenuation"])
			theTrapMutationHandler->setAttenuation(false);
		if (!map["heterogeneous"])
			theNucleationHandler = std::make_shared<
					core::modified::DummyNucleationHandler>();

		// Get the number of dimensions
		int dim = opts.getDimensionNumber();

		// Set-up the grain boundaries from the options
		std::string gbString = opts.getGbString();
		// Build an input stream from the GB string.
		util::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared<std::istringstream>(gbString);
		reader.setInputStream(argSS);
		// Break the string into tokens.
		auto tokens = reader.loadLine();
		// Loop on them
		for (int i = 0; i < tokens.size(); i++) {
			// Switch on the type of grain boundaries
			if (tokens[i] == "X") {
				auto GBAdvecHandler = std::make_shared<
						core::advection::XGBAdvectionHandler>();
				GBAdvecHandler->setLocation(
						strtod(tokens[i + 1].c_str(), NULL));
				GBAdvecHandler->setDimension(dim);
				theAdvectionHandler.push_back(GBAdvecHandler);
			} else if (tokens[i] == "Y") {
				if (dim < 2)
					// A Y grain boundary cannot be used in 1D.
					throw std::string(
							"\nA Y grain boundary CANNOT be used in 1D. Switch to 2D or 3D or remove it.");

				auto GBAdvecHandler = std::make_shared<
						core::advection::YGBAdvectionHandler>();
				GBAdvecHandler->setLocation(
						strtod(tokens[i + 1].c_str(), NULL));
				GBAdvecHandler->setDimension(dim);
				theAdvectionHandler.push_back(GBAdvecHandler);
			} else if (tokens[i] == "Z") {
				if (dim < 3)
					// A Z grain boundary cannot be used in 1D/2D.
					throw std::string(
							"\nA Z grain boundary CANNOT be used in 1D/2D. Switch to 3D or remove it.");

				auto GBAdvecHandler = std::make_shared<
						core::advection::ZGBAdvectionHandler>();
				GBAdvecHandler->setLocation(
						strtod(tokens[i + 1].c_str(), NULL));
				GBAdvecHandler->setDimension(dim);
				theAdvectionHandler.push_back(GBAdvecHandler);
			} else {
				// Wrong GB type
				throw std::string(
						"\nThe type of grain boundary is not known: \""
								+ tokens[i] + "\"");
			}

			i++;
		}

		return;
	}

	/**
	 * Return the flux handler.
	 *
	 *  @return The flux handler.
	 */
	std::shared_ptr<core::flux::IFluxHandler> getFluxHandler() const {
		return theFluxHandler;
	}

	/**
	 * Return the advection handler.
	 *
	 *  @return The advection handler.
	 */
	std::vector<std::shared_ptr<core::advection::IAdvectionHandler> > getAdvectionHandler() const {
		return theAdvectionHandler;
	}

	/**
	 * Return the diffusion handler.
	 *
	 *  @return The diffusion handler.
	 */
	std::shared_ptr<core::diffusion::IDiffusionHandler> getDiffusionHandler() const {
		return theDiffusionHandler;
	}

	/**
	 * Return the modified trap-mutation handler.
	 *
	 *  @return The trap-mutation handler.
	 */
	std::shared_ptr<core::modified::ITrapMutationHandler> getTrapMutationHandler() const {
		return theTrapMutationHandler;
	}

	/**
	 * Return the heterogeneous nucleation handler.
	 *
	 *  @return The nucleation handler.
	 */
	std::shared_ptr<core::modified::IHeterogeneousNucleationHandler> getNucleationHandler() const {
		return theNucleationHandler;
	}
};

} // end namespace material
} // end namespace factory
} // end namespace xolotl

#endif // MATERIALHANDLERFACTORY_H
