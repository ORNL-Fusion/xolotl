// Includes
#include "CvsXYDataProvider.h"
#include <algorithm>

using namespace xolotlViz;

CvsXYDataProvider::CvsXYDataProvider(std::string name) : DataProvider(name) {
}

CvsXYDataProvider::~CvsXYDataProvider() {
}

std::vector<double> CvsXYDataProvider::getAxis1Vector() const {
	std::vector<double> xVector;

	// Loop on all the points in the data vector
	for (auto it = dataPoints->begin();
			it != dataPoints->end(); it++) {

		// Fill the xVector
		addValue(xVector, (*it).x);
	}

	// Add the last point for the mesh (data are in cells, we need one more mesh points than cell points)
	xVector.push_back(xVector.size());

	return xVector;
}

std::vector<double> CvsXYDataProvider::getAxis2Vector() const {
	std::vector<double> yVector;

	// Loop on all the points in the data vector
	for (auto it = dataPoints->begin();
			it != dataPoints->end(); it++) {

		// Fill the yVector
		addValue(yVector, (*it).y);
	}

	// Add the last point for the mesh (data are in cells, we need one more mesh points than cell points)
	yVector.push_back(yVector.size());

	return yVector;
}

std::vector<double> CvsXYDataProvider::getAxis3Vector() const {
	std::vector<double> concentrationVector;

	// Loop on all the points in the data vector
	for (auto it = dataPoints->begin();
			it != dataPoints->end(); it++) {

		// Fill the concentrationVector
		concentrationVector.push_back(std::max((*it).value, 1.0e-16));
	}

	return concentrationVector;
}

void CvsXYDataProvider::addValue(std::vector<double>& vector, double value) const {
	// Check if the value is already in the vector
	auto it = std::find (vector.begin(), vector.end(), value);
	// If it is not, add the value to the vector
	if (it == vector.end()) vector.push_back(value);

	return;
}
