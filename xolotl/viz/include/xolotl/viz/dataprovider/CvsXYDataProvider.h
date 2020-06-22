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
private:
	/**
	 * Declare the constructor as private to force the use of a name
	 */
	CvsXYDataProvider() : DataProvider()
	{
	}

public:
	/**
	 * The default constructor
	 */
	CvsXYDataProvider(const std::string& name);

	/**
	 * The destructor
	 */
	~CvsXYDataProvider();

	/**
	 * Method returning a vector containing the 'x' field of the collection of
	 * DataPoint of the DataProvider.
	 * @return The vector of DataPoint x value.
	 */
	virtual std::vector<double>
	getAxis1Vector() const;

	/**
	 * Method returning a vector containing the 'y' field of the collection of
	 * DataPoint of the DataProvider.
	 * @return The vector of DataPoint y value.
	 */
	virtual std::vector<double>
	getAxis2Vector() const;

	/**
	 * Method returning a vector containing the 'Value' field of the collection
	 * of DataPoint of the DataProvider.
	 * @return The vector of DataPoint value.
	 */
	virtual std::vector<double>
	getAxis3Vector() const;

	/**
	 * Method adding value to a vector if it is not already in the vector.
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
