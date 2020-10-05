#include <fstream>
#include <iostream>
#include <stdexcept>

#include <xolotl/core/Constants.h>
#include <xolotl/core/temperature/HeatEquation1DHandler.h>
#include <xolotl/core/temperature/HeatEquation2DHandler.h>
#include <xolotl/core/temperature/HeatEquation3DHandler.h>
#include <xolotl/core/temperature/TemperatureConstantHandler.h>
#include <xolotl/core/temperature/TemperatureGradientHandler.h>
#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/core/temperature/TemperatureProfileHandler.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
namespace detail
{
auto defaultTemperatureHandlerRegistration =
	xolotl::factory::temperature::TemperatureHandlerFactory::Registration<
		TemperatureHandler>("");
}

TemperatureHandler::TemperatureHandler(const options::Options& options)
{
	// Get the current process ID
	int procId;
	auto xolotlComm = util::getMPIComm();
	MPI_Comm_rank(xolotlComm, &procId);

	if (options.useTemperatureProfileHandlers()) {
		auto tempFileName = options.getTempProfileFilename();
		_strategy = std::make_shared<TemperatureProfileHandler>(tempFileName);
		if (procId == 0)
			std::cout
				<< "TemperatureHandler: Using the time profile defined in: "
				<< tempFileName << std::endl;
	}
	else if (options.useHeatEquationHandlers()) {
		if (util::equal(options.getConstTemperature(), 0.0)) {
			// We are to use a constant temperature handler because the flux is
			// 0.0
			_strategy = std::make_shared<TemperatureConstantHandler>(
				options.getBulkTemperature());
			if (procId == 0)
				std::cout
					<< "TemperatureHandler: Using a constant temperature of: "
					<< options.getBulkTemperature() << " K" << std::endl;
		}
		else {
			// Actually using the heat equation
			// Switch on the dimension
			switch (options.getDimensionNumber()) {
			case 1:
				_strategy = std::make_shared<HeatEquation1DHandler>(
					options.getConstTemperature(),
					options.getBulkTemperature());
				break;
			case 2:
				_strategy = std::make_shared<HeatEquation2DHandler>(
					options.getConstTemperature(),
					options.getBulkTemperature());
				break;
			case 3:
				_strategy = std::make_shared<HeatEquation3DHandler>(
					options.getConstTemperature(),
					options.getBulkTemperature());
				break;
			default:
				// The asked dimension is not good (e.g. -1, 4)
				throw std::runtime_error(
					"\nxolotlFactory: Bad dimension for the heat "
					"equation handler.");
			}

			// Set the heat coefficient which depends on the material
			auto problemType = options.getMaterial();
			// PSI case
			if (problemType == "W100" || problemType == "W110" ||
				problemType == "W111" || problemType == "W211" ||
				problemType == "TRIDYN") {
				_strategy->setHeatCoefficient(core::tungstenHeatCoefficient);
				_strategy->setHeatConductivity(core::tungstenHeatConductivity);
			}
			// NE case
			else if (problemType == "Fuel") {
				_strategy->setHeatCoefficient(core::uo2HeatCoefficient);
				_strategy->setHeatConductivity(core::uo2HeatConductivity);
			}
			// Fe case
			else if (problemType == "Fe") {
				_strategy->setHeatCoefficient(core::feHeatCoefficient);
				_strategy->setHeatConductivity(core::feHeatConductivity);
			}
			if (procId == 0)
				std::cout << "TemperatureHandler: Using the heat equation with "
							 "a flux of: "
						  << options.getConstTemperature()
						  << " W nm-2, and a bulk temperature of: "
						  << options.getBulkTemperature() << " K" << std::endl;
		}
	}
	else if (options.useGradientTemperatureHandlers()) {
		auto temp = options.getConstTemperature();
		// Check if we want a temperature gradient
		double bulkTemp = options.getBulkTemperature();
		// Use a temperature gradient
		_strategy =
			std::make_shared<TemperatureGradientHandler>(temp, bulkTemp);
		if (procId == 0)
			std::cout << "TemperatureHandler: Using a temperature gradient "
						 "with a surface temperature of: "
					  << options.getConstTemperature()
					  << " K, and a bulk temperature of: "
					  << options.getBulkTemperature() << " K" << std::endl;
	}
	else {
		auto temp = options.getConstTemperature();
		// we are to use a constant temperature handler
		_strategy = std::make_shared<TemperatureConstantHandler>(temp);
		if (procId == 0)
			std::cout << "TemperatureHandler: Using a constant temperature of: "
					  << options.getConstTemperature() << " K" << std::endl;
	}
}
} // namespace temperature
} // namespace core
} // namespace xolotl
