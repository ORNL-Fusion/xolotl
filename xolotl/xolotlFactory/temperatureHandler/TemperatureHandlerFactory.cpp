#include <TemperatureHandlerFactory.h>
#include <TemperatureHandler.h>
#include <TemperatureProfileHandler.h>
#include <TemperatureGradientHandler.h>
#include <HeatEquationHandler.h>
#include <MathUtils.h>
#include <fstream>
#include <iostream>
#include <mpi.h>

namespace xolotlFactory {

static std::shared_ptr<xolotlCore::ITemperatureHandler> theTemperatureHandler;

// Create the desired type of handler registry.
bool initializeTempHandler(const xolotlCore::Options &options) {
	// Get the current process ID
	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);

	bool ret = true;

	if (options.useConstTemperatureHandlers()
			&& options.useTemperatureProfileHandlers()) {
		// Only print the error message once when running in parallel
		if (procId == 0) {
			// A constant temperature value AND a temperature profile cannot both be given.
			throw std::string(
					"\nA constant temperature value AND a temperature file cannot both be given.");
		}
	} else if (options.useConstTemperatureHandlers()) {
		auto temp = options.getConstTemperature();
		// Check if we want a temperature gradient
		double bulkTemp = options.getBulkTemperature();
		if (xolotlCore::equal(bulkTemp, 0.0)) {
			// we are to use a constant temperature handler
			theTemperatureHandler = std::make_shared<
					xolotlCore::TemperatureHandler>(temp);
		} else {
			// Use a temperature gradient
			theTemperatureHandler = std::make_shared<
					xolotlCore::TemperatureGradientHandler>(temp, bulkTemp);
		}
	} else if (options.useTemperatureProfileHandlers()) {
		auto tempFileName = options.getTempProfileFilename();
		theTemperatureHandler = std::make_shared<
				xolotlCore::TemperatureProfileHandler>(tempFileName);
	} else if (options.useHeatEquationHandlers()) {
		if (xolotlCore::equal(options.getConstTemperature(), 0.0)) {
			// We are to use a constant temperature handler because the flux is 0.0
			theTemperatureHandler = std::make_shared<
					xolotlCore::TemperatureHandler>(options.getBulkTemperature());
		} else {
			// Actually using the heat equation
			theTemperatureHandler = std::make_shared<
					xolotlCore::HeatEquationHandler>(
					options.getConstTemperature(),
					options.getBulkTemperature());

			// Set the heat coefficient which depends on the material
			auto problemType = options.getMaterial();
			// PSI case
			if (problemType == "W100" || problemType == "W110"
					|| problemType == "W111" || problemType == "W211"
					|| problemType == "TRIDYN") {
				theTemperatureHandler->setHeatCoefficient(
						xolotlCore::tungstenHeatCoefficient);
				theTemperatureHandler->setHeatConductivity(
						xolotlCore::tungstenHeatConductivity);
			}
			// NE case
			else if (problemType == "Fuel") {
				theTemperatureHandler->setHeatCoefficient(
						xolotlCore::uo2HeatCoefficient);
				theTemperatureHandler->setHeatConductivity(
						xolotlCore::uo2HeatConductivity);
			}
			// Fe case
			else if (problemType == "Fe") {
				theTemperatureHandler->setHeatCoefficient(
						xolotlCore::feHeatCoefficient);
				theTemperatureHandler->setHeatConductivity(
						xolotlCore::feHeatConductivity);
			}
		}
	} else {
		// Only print the error message once when running in parallel
		if (procId == 0) {
			std::cerr
					<< "\nWarning: Temperature information has not been given.  Defaulting to constant"
							" temperature = 1000K " << std::endl;
		}
		auto temp = options.getConstTemperature();
		// we are to use a constant temperature handler
		theTemperatureHandler =
				std::make_shared<xolotlCore::TemperatureHandler>(temp);
	}

	return ret;
}

// Provide access to our handler registry.
std::shared_ptr<xolotlCore::ITemperatureHandler> getTemperatureHandler() {
	if (!theTemperatureHandler) {
		// We have not yet been initialized.
		throw std::string("\nxolotlFactory temperature handler requested but "
				"it has not been initialized.");
	}
	return theTemperatureHandler;
}

} // end namespace xolotlFactory

