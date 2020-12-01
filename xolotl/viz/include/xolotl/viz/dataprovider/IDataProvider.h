#ifndef IDATAPROVIDER_H
#define IDATAPROVIDER_H

// Includes
#include <memory>
#include <string>
#include <vector>

#include <xolotl/viz/dataprovider/DataPoint.h>

namespace xolotl
{
namespace viz
{
namespace dataprovider
{
/**
 * IDataProvider describes the structure needed as a link between Xololt outputs
 * and quantities that the user wants to plot.
 */
class IDataProvider
{
public:
	/**
	 * The destructor
	 */
	virtual ~IDataProvider()
	{
	}

	/**
	 * Returns a collection of the data points.
	 *
	 * @return The collection of DataPoint.
	 */
	virtual std::shared_ptr<std::vector<DataPoint>>
	getDataPoints() const = 0;

	/**
	 * Returns the value of the mean of all the data points.
	 *
	 * @return The mean of all the values of the data.
	 */
	virtual double
	getDataMean() const = 0;

	/**
	 * Method filling the data collection.
	 *
	 * @param points A collection of DataPoint.
	 */
	virtual void
	setDataPoints(std::shared_ptr<std::vector<DataPoint>> points) = 0;

	/**
	 * Method returning the vector of quantity that will be plotted on the X
	 * axis. Quantity being x, y, z, t, or value.
	 *
	 * @return The vector of quantity.
	 */
	virtual std::vector<double>
	getAxis1Vector() const = 0;

	/**
	 * Method returning the vector of quantity that will be plotted on the Y
	 * axis. Quantity being x, y, z, t, or value.
	 *
	 * @return The vector of quantity.
	 */
	virtual std::vector<double>
	getAxis2Vector() const = 0;

	/**
	 * Method returning the vector of quantity that will be plotted on the Z
	 * axis. Quantity being x, y, z, t, or value.
	 *
	 * @return The vector of quantity.
	 */
	virtual std::vector<double>
	getAxis3Vector() const = 0;

	/**
	 * Method returning the vector of time steps that will be used for the
	 * VideoPlots.
	 *
	 * @return The vector of time steps.
	 */
	virtual std::vector<double>
	getAxis4Vector() const = 0;

	/**
	 * This operation returns the name of the data.
	 *
	 * @return the name
	 */
	virtual std::string
	getDataName() const = 0;

	/**
	 * This operation sets the name of the data.
	 *
	 * @param name the name of the data
	 */
	virtual void
	setDataName(const std::string& name) = 0;
};

// end class IDataProvider

} /* namespace dataprovider */
} /* namespace viz */
} /* namespace xolotl */
#endif
