#ifndef XCORE_HDF5OBJECT_H
#define XCORE_HDF5OBJECT_H

#include <string>
#include "hdf5.h"

namespace xolotlCore {

// Base class for HDF5 files and objects within them.
class HDF5Object {
private:
    // Name of the object.
    std::string name;

    // HDF5 id associated with the object.
    hid_t id;

protected:
    /**
     * Set the id associated with the object.
     *
     * @param _id The id to be associated with the object.
     */
    void setId(hid_t _id)   { id = _id; }

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
    HDF5Object(std::string _name,
                hid_t _id = H5I_INVALID_HID)
      : name(_name),
        id(_id)
    { }

    ~HDF5Object(void) {
        // Assumes that a derived class closes the id in whatever way
        // is appropriate.
        id = H5I_INVALID_HID;
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


// An HDF5Object that retains the location used when constructing it.
class LocatedHDF5Object : public HDF5Object {
private:
    // Location of the object.
    // Think of it as the "parent" of the object.  In quotes because it 
    // might not be the actual next level up in the namespace within 
    // the HDF5 object.  For example, it might represent the file,
    // and the 'name' would indicate a multi-component path within
    // the file's namespace.
    const HDF5Object& location;

public:
    /**
     * Create a LocatedHDF5Object.
     * Default and copy constructors explicitly disallowed.
     */
    LocatedHDF5Object(void) = delete;
    LocatedHDF5Object(const LocatedHDF5Object& other) = delete;

    /**
     * Create a LocatedHDF5Object.
     *
     * @param _location The location of the object.
     * @param _name Name of the object.
     * @param _id The HDF5 id associated with the object.
     */
    LocatedHDF5Object(const HDF5Object& _location,
                std::string _name,
                hid_t _id = H5I_INVALID_HID)
      : HDF5Object(_name, _id),
        location(_location)
    { }


    /**
     * Obtain the location of the object.
     *
     * @return The location used when defining the object.
     */
    const HDF5Object& getLocation(void) const { return location; }
};

} /* namespace xolotlCore */

#endif // XCORE_HDF5OBJECT_H

