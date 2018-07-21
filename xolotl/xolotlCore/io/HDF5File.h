#ifndef XCORE_HDF5FILE_H
#define XCORE_HDF5FILE_H

#include <string>
#include <vector>
#include "xolotlCore/io/Filesystem.h"
#include "xolotlCore/io/HDF5Object.h"


namespace xolotlCore {

class HDF5File : public HDF5Object {
public:
    // File create/open access modes.
    enum class AccessMode
    {
        OpenReadOnly,
        OpenReadWrite,
        CreateOrTruncateIfExists,
        CreateOrFailIfExists
    };

    // An HDF5 property list.
    class PropertyList : public HDF5Object {
    public:
        /**
         * Create a Property list.
         * Default and copy constructors explicitly disallowed.
         *
         * @param classId HDF5 class ID associated with the property list.
         */
        PropertyList(void) = delete;
        PropertyList(const PropertyList& other) = delete;
        PropertyList(hid_t cls_id)
          : HDF5Object("PropertyList")
        {
            id = H5Pcreate(cls_id);
        }

        /**
         * Release an open property list.
         */
        ~PropertyList(void)
        {
            H5Pclose(id);
        }
    };


    // A group in the HDF5 file.
    class Group : public HDF5Object {
    private:
        // "Location" of our group.
        // Often, it is our parent object, but needn't be if our path
        // has more than one component.
        const HDF5Object& location;

    public:                
        /**
         * Create or open a group within the file.
         * Default and copy constructors explicitly disallowed.
         *
         * @param _location Location within HDF5 file's "filesystem"
         * @param _path Path from _location of group to create/open.
         * @param _create Whether to create group.
         */
        Group(void) = delete;
        Group(const Group& other) = delete;
        Group(const HDF5Object& _location,
                fs::path _path,
                bool _create = true);

        /**
         * Release an open group.
         */
        ~Group(void)
        {
            H5Gclose(id);
        }

        /**
         * Access "location" we were created/opened from.
         *
         * @return The HDF5 object used as "location" when we were 
         * created/opened.
         */
        const HDF5Object& getLocation(void) const   { return location; }
    };


private:
    /**
     * Determine HDF5 flag value that corresponds to given access mode.
     *
     * @param mode An AccessMode constant.
     * @return HDF5 flag that corresponds to the mode parameter.
     */
    static unsigned int toHDF5AccessMode(AccessMode mode);


	/**
	 * Callback for H5Ewalk that constructs an error string for a 
     * specific level of the HDF5 context.
 	 *
	 * @param n Level number within the HDF5 error context.
     * @param eptr The HDF5 error context for level n.
     * @param client_data Data given by application to H5Ewalk.
     * @return Error code indicating whether string was constructed 
     *          successfully.
     */
	static herr_t BuildCurrStackLevelErrorString(uint32_t n,
											const H5E_error2_t* eptr,
											void* client_data);

    /**
     * Construct an error string using the HDF5 context.
     *
     * @return An error string noting the HDF5 context of the error.
     */
    static std::string BuildHDF5ErrorString(void);


public:
    /**
     * Create or open an HDF5 file.
     *
     * @param path Path of file to create or open.
     * @param mode Access mode for creating/opening the file.
     * @param par Whether to access the file with parallel I/O.
     */
    HDF5File(fs::path path, AccessMode mode, bool par = true );

    /**
     * Close the file if open and destroy the in-memory object.
     */
    virtual ~HDF5File(void);
};

} /* namespace xolotlCore */

#endif // XCORE_HDF5FILE_H

