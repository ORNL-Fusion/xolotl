#ifndef ADVECTIONHANDLER_H
#define ADVECTIONHANDLER_H

// Includes
#include "IAdvectionHandler.h"
#include <Constants.h>

namespace xolotlCore {

/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium clusters. It needs to have subclasses
 * that implement the initialize method.
 */
class AdvectionHandler: public IAdvectionHandler {
protected:

	//! The location of the sink
	double location;

	//! The vector containing the indices of the advecting clusters
	std::vector<int> indexVector;

	//! The vector containing the value of the sink strength (called A) of the advecting clusters
	std::vector<double> sinkStrengthVector;

	//! The number of dimensions of the problem
	int dimension;

public:

	//! The Constructor
	AdvectionHandler() : location(0.0), dimension(0) {}

	//! The Destructor
	~AdvectionHandler() {}

	/**
	 * Set the number of dimension.
	 *
	 * @param dim The number of dimensions of the problem
	 */
	void setDimension(int dim) {dimension = dim;}

	/**
	 * Set the location of the sink.
	 *
	 * @param pos The location of the sink
	 */
	void setLocation(double pos) {location = pos;}

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	int getNumberOfAdvecting() {return indexVector.size();}

	/**
	 * Get the vector of index of advecting clusters in the network.
	 *
	 * @return The vector of advecting clusters
	 */
	std::vector<int> getIndexVector() {return indexVector;}

	/**
	 * Get the location of the sink.
	 *
	 * @return The location
	 */
	double getLocation() {return location;}

};
//end class AdvectionHandler

} /* end namespace xolotlCore */
#endif
