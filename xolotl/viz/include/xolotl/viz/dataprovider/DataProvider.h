#ifndef DATAPROVIDER_H
#define DATAPROVIDER_H

// Includes
#include <xolotl/util/Identifiable.h>
#include <xolotl/viz/dataprovider/IDataProvider.h>

namespace xolotl
{
namespace viz
{
namespace dataprovider
{
/**
 *  Realization of the IDataProvider interface. This is a general class with
 * general methods, to actually get data from the data provider, one needs to
 * use the subclasses.
 */
class DataProvider : public virtual IDataProvider, public util::Identifiable
{
protected:
	/**
	 * Declare the constructor as private to force the use of a name
	 */
	DataProvider() : util::Identifiable("unused")
	{
	}

	/**
	 * Collection of data points.
	 */
	std::shared_ptr<std::vector<DataPoint>> dataPoints;

public:
	/**
	 * Name of the data contained in the data attribute.
	 */
	std::string dataName;

	/**
	 * The default constructor.
	 */
	DataProvider(const std::string& name);

	/**
	 * The destructor.
	 */
	~DataProvider();

	/**
	 * Returns a collection of the data points.
	 * \see IDataProvider.h
	 */
	std::shared_ptr<std::vector<DataPoint>>
	getDataPoints() const;

	/**
	 * Method filling the data collection.
	 * \see IDataProvider.h
	 */
	void
	setDataPoints(std::shared_ptr<std::vector<DataPoint>> points);

	/**
	 * Returns the value of the mean of all the data points.
	 * \see IDataProvider.h
	 */
	double
	getDataMean() const;

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of quantity that will be plotted on the X axis.
	 * Quantity being x, y, z, t, or value.
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis1Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of quantity that will be plotted on the Y axis.
	 * Quantity being x, y, z, t, or value.
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis2Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of quantity that will be plotted on the Z axis.
	 * Quantity being x, y, z, t, or value.
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis3Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * Method that has to be overwritten by subclasses.
	 * Should return the vector of time steps that will be used for the
	 * VideoPlots. \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis4Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * This operation returns the name of the data.
	 * @return the name
	 */
	virtual std::string
	getDataName() const
	{
		return dataName;
	}

	/**
	 * This operation sets the name of the data.
	 * @param name the name of the data
	 */
	virtual void
	setDataName(const std::string& name)
	{
		dataName = name;
	}
};

// end class DataProvider

} /* namespace dataprovider */
} /* namespace viz */
} /* namespace xolotl */
#endif
