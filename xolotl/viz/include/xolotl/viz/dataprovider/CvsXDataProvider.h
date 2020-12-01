#ifndef CVSXDATAPROVIDER_H
#define CVSXDATAPROVIDER_H

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
class CvsXDataProvider : public DataProvider
{
public:
	/**
	 * The default constructor
	 */
	CvsXDataProvider() : DataProvider()
	{
	}

	/**
	 * The destructor
	 */
	~CvsXDataProvider()
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
	 * @return The vector of DataPoint value.
	 *
	 * \see IDataProvider.h
	 */
	virtual std::vector<double>
	getAxis2Vector() const;
};

// end class CvsXDataProvider

} /* namespace dataprovider */
} /* namespace viz */
} /* namespace xolotl */

#endif
