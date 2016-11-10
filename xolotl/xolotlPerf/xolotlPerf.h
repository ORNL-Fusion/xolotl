#ifndef XOLOTLPERF_H
#define XOLOTLPERF_H

#include <memory>
#include <sstream>
#include "IHandlerRegistry.h"
#include "RuntimeError.h"

namespace xolotlPerf {

/**
 * Detect the type of performance handlers registry to create based
 * on a string argument (e.g., taken from the command line).
 * Throws an std::invalid_argument exception if arg does not
 * specify a registry type we know about.
 *
 * NOTE: We define the function in the header to avoid adding a dependency
 * on the performance library for the command line library which uses
 * this function to parse its command line.
 *
 * @param arg String description of performance handler registry to create.
 * @return The performance handler registry type to create.
 */
inline IHandlerRegistry::RegistryType toPerfRegistryType(
		const std::string& arg) {
	IHandlerRegistry::RegistryType ret;

	if (arg == "dummy") {
		ret = IHandlerRegistry::dummy;
	} else if (arg == "std") {
		ret = IHandlerRegistry::std;
	} else if (arg == "os") {
		ret = IHandlerRegistry::os;
	} else if (arg == "papi") {
		ret = IHandlerRegistry::papi;
	} else {
		std::ostringstream estr;
		estr << "Invalid performance handler argument \"" << arg << "\" seen.";
		throw std::invalid_argument(estr.str());
	}
	return ret;
}

/**
 * Initialize the performance library for using the desired type of handlers.
 * Throws a std::invalid_argument if caller requests a registry type that
 * we do not support.
 *
 * @param rtype Type of handlerRegistry to create.
 */
void initialize(IHandlerRegistry::RegistryType rtype);

/**
 * Access the handler registry.
 * Throws a std::runtime_error if called before the xolotlPerf classes
 * have been initialized.
 *
 *  @return The handler registry object.
 */
std::shared_ptr<IHandlerRegistry> getHandlerRegistry(void);

} // end namespace xolotlPerf

#endif // XOLOTLPERF_H
