#ifndef MATERIALHANDLERFACTORY_H
#define MATERIALHANDLERFACTORY_H

#include <memory>
#include <IMaterialFactory.h>
#include <DummyDiffusionHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>
#include <DummyReSolutionHandler.h>
#include <DummyNucleationHandler.h>
#include <TokenizedLineReader.h>
#include <XGBAdvectionHandler.h>
#include <YGBAdvectionHandler.h>
#include <ZGBAdvectionHandler.h>

namespace xolotlFactory {

/**
 * Realizes the IMaterialFactory interface. Handles the flux, the advection, modified
 * trap mutation, and for a specific material.
 */
class MaterialFactory: public IMaterialFactory {
protected:

	//! The flux handler
	std::shared_ptr<xolotlCore::IFluxHandler> theFluxHandler;

	//! The advection handler
	std::vector<std::shared_ptr<xolotlCore::IAdvectionHandler> > theAdvectionHandler;

	//! The diffusion handler
	std::shared_ptr<xolotlCore::IDiffusionHandler> theDiffusionHandler;

	//! The modified trap-mutation handler
	std::shared_ptr<xolotlCore::ITrapMutationHandler> theTrapMutationHandler;

	//! The re-solution handler
	std::shared_ptr<xolotlCore::IReSolutionHandler> theReSolutionHandler;

	//! The heterogeneous nucleation handler
	std::shared_ptr<xolotlCore::IHeterogeneousNucleationHandler> theNucleationHandler;

public:

	/**
	 * The constructor creates the handlers.
	 */
	MaterialFactory() {
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
	virtual void initializeMaterial(const xolotlCore::Options &options) {
		// Wrong if both he flux and time profile options are used
		if (options.useFluxAmplitude() && options.useFluxTimeProfile()) {
			// A constant flux value AND a time profile cannot both be given.
			throw std::string(
					"\nA constant flux value AND a time profile cannot both be given.");
		} else if (options.useFluxAmplitude()) {
			// Set the constant value of the flux
			theFluxHandler->setFluxAmplitude(options.getFluxAmplitude());
		} else if (options.useFluxTimeProfile()) {
			// Initialize the time profile
			theFluxHandler->initializeTimeProfile(options.getFluxProfileName());
		}

		// Get the process map
		auto map = options.getProcesses();
		// Set dummy handlers when needed
		if (!map["diff"])
			theDiffusionHandler = std::make_shared<
					xolotlCore::DummyDiffusionHandler>();
		if (!map["advec"]) {
			// Clear the advection handler
			theAdvectionHandler.clear();
			// To replace it by a dummy one
			theAdvectionHandler.push_back(
					std::make_shared<xolotlCore::DummyAdvectionHandler>());
		}
		if (!map["modifiedTM"])
			theTrapMutationHandler = std::make_shared<
					xolotlCore::DummyTrapMutationHandler>();
		if (!map["attenuation"])
			theTrapMutationHandler->setAttenuation(false);
		if (!map["resolution"])
			theReSolutionHandler = std::make_shared<
					xolotlCore::DummyReSolutionHandler>();
		if (!map["heterogeneous"])
			theNucleationHandler = std::make_shared<
					xolotlCore::DummyNucleationHandler>();

		// Get the number of dimensions
		int dim = options.getDimensionNumber();

		// Set-up the grain boundaries from the options
		std::string gbString = options.getGbString();
		// Build an input stream from the GB string.
		xolotlCore::TokenizedLineReader<std::string> reader;
		auto argSS = std::make_shared<std::istringstream>(gbString);
		reader.setInputStream(argSS);
		// Break the string into tokens.
		auto tokens = reader.loadLine();
		// Loop on them
		for (int i = 0; i < tokens.size(); i++) {
			// Switch on the type of grain boundaries
			if (tokens[i] == "X") {
				auto GBAdvecHandler = std::make_shared<
						xolotlCore::XGBAdvectionHandler>();
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
						xolotlCore::YGBAdvectionHandler>();
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
						xolotlCore::ZGBAdvectionHandler>();
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
	std::shared_ptr<xolotlCore::IFluxHandler> getFluxHandler() const {
		return theFluxHandler;
	}

	/**
	 * Return the advection handler.
	 *
	 *  @return The advection handler.
	 */
	std::vector<std::shared_ptr<xolotlCore::IAdvectionHandler> > getAdvectionHandler() const {
		return theAdvectionHandler;
	}

	/**
	 * Return the diffusion handler.
	 *
	 *  @return The diffusion handler.
	 */
	std::shared_ptr<xolotlCore::IDiffusionHandler> getDiffusionHandler() const {
		return theDiffusionHandler;
	}

	/**
	 * Return the modified trap-mutation handler.
	 *
	 *  @return The trap-mutation handler.
	 */
	std::shared_ptr<xolotlCore::ITrapMutationHandler> getTrapMutationHandler() const {
		return theTrapMutationHandler;
	}

	/**
	 * Return the Re-solution handler.
	 *
	 *  @return The re-solution handler.
	 */
	std::shared_ptr<xolotlCore::IReSolutionHandler> getReSolutionHandler() const {
		return theReSolutionHandler;
	}

	/**
	 * Return the heterogeneous nucleation handler.
	 *
	 *  @return The nucleation handler.
	 */
	std::shared_ptr<xolotlCore::IHeterogeneousNucleationHandler> getNucleationHandler() const {
		return theNucleationHandler;
	}
};

} // end namespace xolotlFactory

#endif // MATERIALHANDLERFACTORY_H
