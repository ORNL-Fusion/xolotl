#ifndef XCORE_HDF5OBJECT_H
#define XCORE_HDF5OBJECT_H

#include <string>
#include "hdf5.h"

namespace xolotlCore {

// Base class for HDF5 files and objects within them.
class HDF5Object {
protected:

    // Name of the object.
    std::string name;

    // HDF5 id associated with the object.
    hid_t id;

public:

    /**
     * Create an HDF5Object.
     * Default and copy constructors explicitly disallowed.
     */
    HDF5Object(void) = delete;
    HDF5Object(const HDF5Object& other) = delete;

    /**
     * Create an HDF5Object.
     *
     * @param _name Name of the object.
     * @param _id The HDF5 id associated with the object.
     */
    HDF5Object(std::string _name, hid_t _id = H5I_INVALID_HID)
      : name(_name),
        id(_id) {

    }

    /**
     * Access the HDF5 id associated with the object.
     *
     * @return The HDF5 id associated with the object.
     */
    hid_t getId(void) const { return id; }


    /**
     * Obtain the name associated with the object.
     *
     * @return The name associated with the object.
     */
    std::string getName(void) const { return name; }
};

} /* namespace xolotlCore */

#endif // XCORE_HDF5OBJECT_H

