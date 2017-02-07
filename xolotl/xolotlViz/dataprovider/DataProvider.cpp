// Includes
#include "DataProvider.h"

using namespace xolotlViz;

DataProvider::DataProvider(const std::string& name) : xolotlCore::Identifiable(name) {
}

DataProvider::~DataProvider() {
}

std::shared_ptr< std::vector<Point> > DataProvider::getDataPoints() const {
	return dataPoints;
}

void DataProvider::setPoints(std::shared_ptr< std::vector<Point> > points) {
	dataPoints = points;
	return;
}

double DataProvider::getDataMean() const {
	// The size of the data vector
	int size = dataPoints->size();

	// Use to add the value of each Point
	double valueSum = 0.0;

	// Loop on all the points in the data vector
	for (auto it = dataPoints->begin();
			it != dataPoints->end(); ++it) {

		// Add the current value to the sum
		valueSum += (*it).value;
	}

	// Result equals the sum divided by the size
	double result = (double) valueSum / size;

	return result;
}
