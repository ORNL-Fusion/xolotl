#ifndef XPERF_RUNTIME_ERROR_H
#define XPERF_RUNTIME_ERROR_H

#include <stdexcept>

namespace xolotl {
namespace perf {

/// A std::runtime_error that also maintains an error code.
/// Exceptions of this type are thrown by some of the xolotlPerf
/// classes when an error code is available, such as when
/// an error occurs in a PAPI function.
class runtime_error: public std::runtime_error {
private:
	/// The error code.
	int code;

public:
	/// Construct a runtime_error.
	///
	/// @param _what The context or error string associated with the exception.
	/// @param _code The error code associated with the exception.
	runtime_error(const std::string& _what, int _code) :
			std::runtime_error(_what), code(_code) {
	}

	/// Retreive the error code.
	///
	/// @return The error code associated with the exception.
	int getCode(void) const {
		return code;
	}
};

}//end namespace perf
}//end namespace xolotl

#endif // XPERF_RUNTIME_ERROR_H
