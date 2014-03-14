#ifndef GPTLHARDWARECOUNTER_H
#define GPTLHARDWARECOUNTER_H

#include <vector>
#include <unordered_map>
#include "IHardwareCounter.h"
#include "HardwareQuantities.h"
#include "Identifiable.h"


namespace xolotlPerf{

struct hardwareQuantity {
	std::string hardwareQuantityString;
	int papiQuantity;
	std::string papiQuantityString;
} ;

/**
 * The GPTLHardwareCounter class is instantiated by the StandardHandlerRegistry class
 * and realizes the IHardwareCounter interface to gather hardware performance counter data
 * found by utilizing the PAPI (Performance Application Programming Interface) library via the
 * General Purpose Timing Library (GPTL).
 */
class GPTLHardwareCounter: public IHardwareCounter, public xolotlCore::Identifiable
{

private:

		/**
		 * The hardware quantities this GPTLHardwareCounter monitors.
		 */
		std::vector<HardwareQuantities> quantities;

		/**
		 * A map that is used to map the hardware quantities to their
		 * corresponding PAPI counterparts
		 */
		static std::unordered_map<int, hardwareQuantity> AllHardwareQuantitiesMap;

		/**
		 * A map that is used to map the hardware quantities to their
		 * corresponding PAPI counterparts
		 */
		std::unordered_map<int, hardwareQuantity> hardwareQuantitiesMap;

		/**
		 * The default constructor is private because HardwareCounters must
		 * always be given a name and a vector of quantities to
		 * be monitored.
		 */
		GPTLHardwareCounter()
			: xolotlCore::Identifiable("unused"),
			  quantities(std::vector<HardwareQuantities>(0))
		{ }


        /**
         * Initialize the map of hardware quantities information.
         */
        static void InitHardwareQuantitiesMap(void);

public:

	/**
	 * IHardwareCounter constructor that takes the name and a
	 * list of the different quantities it should monitor.
	 *
	 * @param counterName The IHardwareCounter's name
	 * @param counterQuantities The vector of quantities the GPTLHardwareCounter
	 * will monitor
	 */
	GPTLHardwareCounter(std::string name,
			const std::vector<HardwareQuantities> &counterQuantities);

	/**
	 * The destructor
	 */
    virtual ~GPTLHardwareCounter() { }

    /**
     * This operation returns a list of values of the, initially specified,
     * hardware quantities monitored by the GPTLHardwareCounter.
     */
    virtual std::vector<double> getValues() const;

	/**
	 * This operation returns the list of hardware
	 * quantities monitored by the GPTLHardwareCounter.
	 */
	virtual std::vector<std::string> getHardwareQuantities() const;

};  //end class GPTLHardwareCounter

}  //end namespace xolotlPerf

#endif
