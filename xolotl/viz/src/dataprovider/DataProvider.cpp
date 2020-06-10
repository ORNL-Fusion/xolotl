// Includes
#include <xolotl/viz/dataprovider/DataProvider.h>

namespace xolotl {
namespace viz {
namespace dataprovider {

DataProvider::DataProvider(const std::string& name) : util::Identifiable(name) {
}

DataProvider::~DataProvider() {
}

std::shared_ptr< std::vector<DataPoint> > DataProvider::getDataPoints() const {
	return dataPoints;
}

void DataProvider::setDataPoints(std::shared_ptr< std::vector<DataPoint> > points) {
	dataPoints = points;
	return;
}

double DataProvider::getDataMean() const {
	// The size of the data vector
	int size = dataPoints->size();

	// Use to add the value of each DataPoint
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

}
}
}
