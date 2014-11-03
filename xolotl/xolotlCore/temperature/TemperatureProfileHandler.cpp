#include "TemperatureProfileHandler.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cmath>

namespace xolotlCore {

void TemperatureProfileHandler::initializeTemperature() {

	// Open file dataFile.dat containing the time and temperature
	std::ifstream inputFile(tempFile.c_str());
	std::string line;

	while (getline(inputFile, line)) {
		if (!line.length() || line[0] == '#')
			continue;
		double xtemp = 0.0, ytemp = 0.0;
		sscanf(line.c_str(), "%lf %lf", &xtemp, &ytemp);
		time.push_back(xtemp);
		temp.push_back(ytemp);
	}

}

double TemperatureProfileHandler::getTemperature(std::vector<double> position,
		double currentTime) const {

	double f;

	// if x is smaller than or equal to xi[0]
	if (currentTime <= time[0])
		return f = temp[0];

	// if x is greater than or equal to xi[n-1]
	if (currentTime >= time[time.size() - 1])
		return f = temp[time.size() - 1];

	// loop to determine the interval x falls in, ie x[k] < x < x[k+1]
	for (int k = 0; k < time.size() - 1; k++) {
		if (currentTime < time[k]) continue;
		if (currentTime > time[k + 1]) continue;

		f = temp[k]
				+ (temp[k + 1] - temp[k]) * (currentTime - time[k])
						/ (time[k + 1] - time[k]);
		break;
	}

	return f;
}

}/* end namespace xolotlCore */

