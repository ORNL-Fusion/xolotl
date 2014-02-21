#ifndef HARDWARECOUNTER_H
#define HARDWARECOUNTER_H

#include <string>
#include <vector>
#include <memory>
#include "HardwareQuantities.h"

using namespace std;

namespace xolotlPerf{

/**
 * Realizations of this interface are responsible for the
 * collection of hardware performance counter data.
 */
class HardwareCounter {

private:

	/**
	 * The default constructor is private because HardwareCounters must
	 * always be initialized with a name and a vector of quantities to
	 * be monitored.
	 */
//	HardwareCounter();

protected:

	/**
	 * The name of this HardwareCounter.
	 */
	std::string name;

	/**
	 * The vector of quantities the HardwareCounter will monitor
	 */
	std::vector<HardwareQuantities> quantities;
	// add MAX_QUANTITIES to enum to initialize the quantities vector to contain
	// MAX_QUANTITIES elements

	/**
	 * The values of this HardwareCounter.
	 */
	std::vector<long long> values;
	// size of values = quantities.size()

public:

	/**
	 * HardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param aname The HardwareCounter's name
	 * @param hquantities The HardwareCounter's list of quantities
	 */
	HardwareCounter(std::string aname, const std::vector<HardwareQuantities> &hquantities);
//    HardwareCounter(std::string name, std::vector<HardwareQuantities> hquantities ) :
//    	name(aname), quantities(hquantities) {}

	/**
	 * The copy constructor.
	 * @param other The HardwareCounter to copy
	 */
	HardwareCounter(const HardwareCounter &other);

	/**
	 * The destructor
	 */
	virtual ~HardwareCounter();

	/**
	 * This operation returns a HardwareCounter that is created using the copy
	 * constructor. If this HardwareCounter is actually a subclass of HardwareCounter, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * HardwareCounter.
	 * @return A copy of this HardwareCounter.
	 */
//	virtual std::shared_ptr<HardwareCounter> clone();

    // This operation returns a list of values of the, initially specified, PAPI preset quantities monitored by the HardwareCounter.
//    virtual std::vector<long long> getValues() const = 0;
    virtual std::vector<long long> getValues() const;

    // This operation returns the name of the HardwareCounter.
    virtual const std::string getName() const;
//    const std::string getName() const {return aname;}

    // This operation increments the HardwareCounter.
    virtual void increment() = 0;

};  //end class HardwareCounter

}  //end namespace xolotlPerf

#endif
