#ifndef W211MATERIALHANDLERFACTORY_H
#define W211MATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <W211FitFluxHandler.h>
#include <W211AdvectionHandler.h>
#include <W211TrapMutationHandler.h>
#include <Diffusion1DHandler.h>
#include <Diffusion2DHandler.h>
#include <Diffusion3DHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for a (211) oriented tungsten material.
 */
class W211MaterialFactory: public MaterialFactory {
private:

	/**
	 * The default constructor is private.
	 */
	W211MaterialFactory() {
	}

public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	W211MaterialFactory(int dim) {
		theFluxHandler = std::make_shared<xolotlCore::W211FitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::W211AdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::W211TrapMutationHandler>();
		theReSolutionHandler = std::make_shared<
				xolotlCore::DummyReSolutionHandler>();

		// Switch on the dimension for the diffusion handler
		switch (dim) {
		case 0:
			theDiffusionHandler = std::make_shared<
					xolotlCore::DummyDiffusionHandler>();
			break;
		case 1:
			theDiffusionHandler = std::make_shared<
					xolotlCore::Diffusion1DHandler>();
			break;
		case 2:
			theDiffusionHandler = std::make_shared<
					xolotlCore::Diffusion2DHandler>();
			break;
		case 3:
			theDiffusionHandler = std::make_shared<
					xolotlCore::Diffusion3DHandler>();
			break;
		default:
			// The asked dimension is not good (e.g. -1, 4)
			throw std::string(
					"\nxolotlFactory: Bad dimension for the W211 material factory.");
		}

		return;
	}

	/**
	 * The destructor
	 */
	~W211MaterialFactory() {
	}

	/**
	 * Initialize the material conditions with the different given options.
	 *
	 * @param options The Xolotl options.
	 */
	void initializeMaterial(const xolotlCore::Options &options) {
		// First call the original one
		MaterialFactory::initializeMaterial(options);

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
};

} // end namespace xolotlFactory

#endif // W211MATERIALHANDLERFACTORY_H
