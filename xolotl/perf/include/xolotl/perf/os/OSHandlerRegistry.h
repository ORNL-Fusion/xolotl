#ifndef OSHANDLERREGISTRY_H
#define OSHANDLERREGISTRY_H

#include <xolotl/perf/standard/StdHandlerRegistry.h>

namespace xolotl {
namespace perf {
namespace os {

/**
 * Factory for building performance data collection objects that 
 * use an OS/runtime timer API.  No hardware performance counters
 * are available.
 */
class OSHandlerRegistry: public standard::StdHandlerRegistry {
public:

	/// Construct a handler registry.
	OSHandlerRegistry(void) {
	}

	/// Destroy the handler registry.
	virtual ~OSHandlerRegistry(void) {
	}

	/**
	 * Look up and return a named timer.
	 * Create the timer if it does not already exist.
	 *
	 * @param name The object's name.
	 * @return The object with the given name.
	 */
	std::shared_ptr<ITimer> getTimer(const std::string& name) override;

	/**
	 * Look up and return a hardware counter set in the current scope.
	 * Create the event counter set if it does not already exist.
	 *
	 * @param name The object's name.
	 * @return The object with the given name.
	 */
	std::shared_ptr<IHardwareCounter> getHardwareCounter(
			const std::string& name,
            const IHardwareCounter::SpecType& ctrSpec) override;
};

} // namespace os
} // namespace perf
} // namespace xolotl

#endif // OSHandlerRegistry
