#ifndef XCORE_HDF5EXCEPTION_H
#define XCORE_HDF5EXCEPTION_H

#include <stdexcept>

namespace xolotlCore
{

class HDF5Exception : public std::runtime_error
{
public:
	/**
	 * Construct an HDF5 exception object.
	 * Default constructor explicitly disallowed.
	 */
    HDF5Exception(void) = delete;

    /**
     * Construct an HDF5 exception with the given message.
     *
     * @param msg The message to associate with the exception.
     */
    HDF5Exception(std::string msg)
      : std::runtime_error(msg)
    { }

    /**
     * Construct an HDF5 exception object as a copy of another.
     *
     * @param other The exception to copy.
     */
    HDF5Exception(const HDF5Exception& other) = default;
};

} // namespace xolotlCore

#endif // XCORE_HDF5EXCEPTION_H
