#ifndef STANDARDHANDLERREGISTRY_H
#define STANDARDHANDLERREGISTRY_H

#include <string>
#include <vector>
#include <memory>
#include "IHandlerRegistry.h"
#include "GPTLTimer.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLTimer
#include "GPTLEventCounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLEventCounter
#include "GPTLHardwareCounter.h" //Dependency Generated Source:StandardHandlerRegistry Target:GPTLHardwareCounter

namespace xolotlPerf {

/*The StandardHandlerRegistry class realizes the interface IHandlerRegistry
 * to acquire the performance data found by implementing the performance interfaces
 * Timer, EventCounter, and HardwareCounter.
 */
class StandardHandlerRegistry : public IHandlerRegistry
{

//    public:
//
//        StandardHandlerRegistry(StandardHandlerRegistry & arg);
////        StandardHandlerRegistry();
//
//        ~StandardHandlerRegistry();
//
//        // This operation returns the Timer specified by the parameter.
//        std::shared_ptr<Timer> getTimer(const std::string name) const;
//
//        // This operation returns the EventCounter specified by the parameter.
//        std::shared_ptr<EventCounter> getEventCounter(const std::string name) const;
//
//        // This operation returns the specified HardwareCounter.
//        std::shared_ptr<HardwareCounter> getHardwareCounter(const std::string name,
//        		const std::shared_ptr<std::vector<HardwareQuantities> > quantities) const;
////        std::shared_ptr<HardwareCounter> getHardwareCounter(const std::string name,
////        			const std::vector<HardwareQuantities> quantities) const;
//
//        // This operation returns a list of values of the, initially specified, PAPI
//        // preset quantities monitored by the HardwareCounter.
//        const std::shared_ptr<std::vector<HardwareQuantities> > getHardwareQuantities() const;
//
//        //        const std::vector<HardwareQuantities> & getHardwareQuantities() const;
//
//        // This operation outputs the information gathered.
//        virtual void dump(std::ostream out) const;

};  //end class StandardHandlerRegistry

}//end namespace xolotlPerf

#endif
