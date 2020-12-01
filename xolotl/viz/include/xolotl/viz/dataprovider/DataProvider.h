#ifndef DATAPROVIDER_H
#define DATAPROVIDER_H

// Includes
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
class DataProvider : public virtual IDataProvider
{
protected:
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
	DataProvider()
	{
	}

	/**
	 * The destructor.
	 */
	~DataProvider()
	{
	}

	/**
	 * \see IDataProvider.h
	 */
	std::shared_ptr<std::vector<DataPoint>>
	getDataPoints() const
	{
		return dataPoints;
	}

	/**
	 * \see IDataProvider.h
	 */
	void
	setDataPoints(std::shared_ptr<std::vector<DataPoint>> points)
	{
		dataPoints = points;
	}

	/**
	 * \see IDataProvider.h
	 */
	double
	getDataMean() const;

	/**
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis1Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis2Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis3Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis4Vector() const
	{
		return std::vector<double>();
	}

	/**
	 * \see IDataProvider.h
	 */
	virtual std::string
	getDataName() const
	{
		return dataName;
	}

	/**
	 * \see IDataProvider.h
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
