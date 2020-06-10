#ifndef CVSXDATAPROVIDER_H
#define CVSXDATAPROVIDER_H

// Includes
#include <vector>
#include <xolotl/viz/dataprovider/DataProvider.h>

namespace xolotl {
namespace viz {
namespace dataprovider {

/**
 * Subclass of DataProvider that will provide the methods to give the value (concentration here)
 * and X data to a ScatterPlot.
 */
class CvsXDataProvider: public DataProvider {

private:

	/**
	 * Declare the constructor as private to force the use of a name
	 */
	CvsXDataProvider() :
			DataProvider() {
	}

public:

	/**
	 * The default constructor
	 */
	CvsXDataProvider(const std::string& name);

	/**
	 * The destructor
	 */
	~CvsXDataProvider();

	/**
	 * Method returning a vector containing the 'x' field of the collection of DataPoint of the DataProvider.
	 * @return The vector of DataPoint x value.
	 */
	virtual std::vector<double> getAxis1Vector() const;

	/**
	 * Method returning a vector containing the 'value' field of the collection of DataPoint of the DataProvider.
	 * @return The vector of DataPoint value.
	 */
	virtual std::vector<double> getAxis2Vector() const;

};

//end class CvsXDataProvider

} /* namespace dataprovider */
} /* namespace viz */
} /* namespace xolotl */

#endif
