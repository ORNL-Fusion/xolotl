#ifndef ADVECTIONHANDLER_H
#define ADVECTIONHANDLER_H

// Includes
#include <xolotl/core/Constants.h>
#include <xolotl/core/advection/IAdvectionHandler.h>

namespace xolotl
{
namespace core
{
namespace advection
{
/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium clusters. It needs to
 * have subclasses that implement the initialize method.
 */
class AdvectionHandler : public IAdvectionHandler
{
protected:
	//! The location of the sink
	double location;

	//! The collection of advecting clusters.
	std::vector<std::size_t> advectingClusters;

	//! The vector containing the value of the sink strength (called A) of the
	//! advecting clusters
	std::vector<double> sinkStrengthVector;

	//! The number of dimensions of the problem
	int dimension;

public:
	//! The Constructor
	AdvectionHandler() : location(0.0), dimension(0)
	{
	}

	//! The Destructor
	~AdvectionHandler()
	{
	}

	/**
	 * Set the number of dimension.
	 *
	 * @param dim The number of dimensions of the problem
	 */
	void
	setDimension(int dim) override
	{
		dimension = dim;
	}

	/**
	 * Set the location of the sink.
	 *
	 * @param pos The location of the sink
	 */
	void
	setLocation(double pos) override
	{
		location = pos;
	}

	/**
	 * Get the total number of advecting clusters in the network.
	 *
	 * @return The number of advecting clusters
	 */
	int
	getNumberOfAdvecting() const override
	{
		return advectingClusters.size();
	}

	/**
	 * Get the vector of index of advecting clusters in the network.
	 *
	 * @return The vector of advecting clusters
	 */
	const std::vector<std::size_t>&
	getAdvectingClusters() override
	{
		return advectingClusters;
	}

	/**
	 * Get the location of the sink.
	 *
	 * @return The location
	 */
	double
	getLocation() const override
	{
		return location;
	}
};
// end class AdvectionHandler

} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
#endif
