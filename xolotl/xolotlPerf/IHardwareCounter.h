#ifndef IHARDWARECOUNTER_H
#define IHARDWARECOUNTER_H

#include "IIdentifiable.h"

using namespace std;

namespace xolotlPerf{

/**
 * Realizations of this interface are responsible for the
 * collection of hardware performance counter data.
 */
class IHardwareCounter : public virtual xolotlCore::IIdentifiable {

public:

	/**
	 * The destructor
	 */
	virtual ~IHardwareCounter() { }

    /**
     * This operation returns a list of values of the, initially specified,
     * hardware quantities monitored by the IHardwareCounter.
     */
    virtual std::vector<double> getValues() const = 0;

	/**
	 * This operation returns the list of hardware
	 * quantities monitored by the IHardwareCounter.
	 */
	virtual std::vector<std::string> getHardwareQuantities() const = 0;

};  //end class IHardwareCounter

}  //end namespace xolotlPerf

#endif
