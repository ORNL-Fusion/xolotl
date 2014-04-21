#include "TemperatureHandler.h"
#include <vector>

using namespace xolotlSolver;

TemperatureHandler::TemperatureHandler(){

}

TemperatureHandler::~TemperatureHandler() {

}

double TemperatureHandler::getTemperature(std::vector<double> position, double currentTime) const {

	// Set the temperature to 1000K
	double temperature = 1000.0;

	return temperature;
}
