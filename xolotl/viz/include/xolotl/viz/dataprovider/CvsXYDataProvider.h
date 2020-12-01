#ifndef CVSXYDATAPROVIDER_H
#define CVSXYDATAPROVIDER_H

// Includes
#include <vector>

#include <xolotl/viz/dataprovider/DataProvider.h>

namespace xolotl
{
namespace viz
{
namespace dataprovider
{
/**
 * Subclass of DataProvider that will provide the methods to give the value
 * (concentration here) and X data to a ScatterPlot.
 */
class CvsXYDataProvider : public DataProvider
{
public:
	/**
	 * The default constructor
	 */
	CvsXYDataProvider() : DataProvider()
	{
	}

	/**
	 * The destructor
	 */
	~CvsXYDataProvider()
	{
	}

	/**
	 * @return The vector of DataPoint x value.
	 *
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis1Vector() const;

	/**
	 * @return The vector of DataPoint y value.
	 *
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis2Vector() const;

	/**
	 * @return The vector of DataPoint value.
	 *
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis3Vector() const;

private:
	/**
	 * Method adding value to a vector if it is not already in the vector.
	 *
	 * @param vector The vector to which to add the value.
	 * @param value The value to add.
	 */
	virtual void
	addValue(std::vector<double>& vector, double value) const;
};

// end class CvsXYDataProvider

} /* namespace dataprovider */
} /* namespace viz */
} /* namespace xolotl */

#endif
