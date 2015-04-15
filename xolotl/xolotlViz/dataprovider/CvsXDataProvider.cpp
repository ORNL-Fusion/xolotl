// Includes
#include "CvsXDataProvider.h"

using namespace xolotlViz;

CvsXDataProvider::CvsXDataProvider(std::string name) : DataProvider(name) {
}

CvsXDataProvider::~CvsXDataProvider() {
}

std::vector<double> CvsXDataProvider::getAxis1Vector() const {
	std::vector<double> xVector;

	// Loop on all the points in the data vector
	for (auto it = dataPoints->begin();
			it != dataPoints->end(); ++it) {

		// Fill the xVector
		xVector.push_back((*it).x);
	}

	return xVector;
}

std::vector<double> CvsXDataProvider::getAxis2Vector() const {
	std::vector<double> concentrationVector;

	// Loop on all the points in the data vector
	for (auto it = dataPoints->begin();
			it != dataPoints->end(); ++it) {

		// Fill the concentrationVector
		concentrationVector.push_back(std::max((*it).value, 1.0e-16));
	}

	return concentrationVector;
}
