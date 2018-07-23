#ifndef XCORE_HDF5FILE_H
#define XCORE_HDF5FILE_H

#include <string>
#include <vector>
#include <array>
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

    // A class to manage lifetime of a dataspace in the underlying file.
    class DataSpace : public HDF5Object
    {
    protected:
        /**
         * Create a dataspace.
         *
         * @param _id The HDF5 id of the dataspace.
         */
        DataSpace(hid_t _id)
          : HDF5Object("DataSpace", _id)
        { }

    public:
        /**
         * Create a dataspace.
         * Default and copy constructor explicitly disallowed.
         */
        DataSpace(void) = delete;
        DataSpace(const DataSpace& other) = delete;

        /**
         * Release an open DataSpace.
         */
        ~DataSpace(void)
        {
            H5Sclose(getId());
        }
    };

    // A scalar DataSpace.
    class ScalarDataSpace : public DataSpace
    {
    public:
        /**
         * Construct a DataSpace for a scalar value.
         */
        ScalarDataSpace(void);
    };

    // A "simple" DataSpace.  I.e., one that can be created using
    // HDF5's H5Ssimple* functions.
    class AttributeBase;
    class DataSetBase;
    template<uint32_t Rank>
    class SimpleDataSpace : public DataSpace
    {
    public:
        // Concise type for dimensions of a SimpleDataSpace.
        using Dimensions = std::array<hsize_t, Rank>;

    private:
        // Dimensions of the dataspace.
        Dimensions dims;

    public:
        SimpleDataSpace(const Dimensions& _dims);
        SimpleDataSpace(const Dimensions& _dims, const Dimensions& _maxDims);
        SimpleDataSpace(const AttributeBase& _attr);
        SimpleDataSpace(const DataSetBase& _dset);

        /**
         * Obtain our dimensions.
         *
         * @return Dimensions of our SimpleDataSpace.
         */
        const Dimensions& getDims(void) const { return dims; }

        /**
         * Change our dimensions to the given dimensions.
         *
         * @param _dims The dimensions we should use.
         */
        void setDims(const Dimensions& _dims);
    };

    class TypeBase : public HDF5Object
    {
    private:
        bool shouldClose;

    public:
        TypeBase(void) = delete;
        TypeBase(std::string _name, hid_t _hid, bool _shouldClose)
          : HDF5Object(_name, _hid),
            shouldClose(_shouldClose)
        { }
        TypeBase(const TypeBase& other) = delete;

        ~TypeBase(void)
        {
            if(shouldClose)
            {
                H5Tclose(id);
            }
            id = H5I_INVALID_HID;
        }
    };

    template<typename T>
    class TypeInFile : public TypeBase
    {
    public:
        TypeInFile(void);
        TypeInFile(const TypeInFile& other) = delete;
    };

    // Partial specialization of TypeInFile for vectors.
    template<typename T>
    class TypeInFile<std::vector<T>> : public TypeBase
    {
    public:
        TypeInFile(void);
        TypeInFile(const TypeInFile& other) = delete;
    };

    template<typename T>
    class TypeInMemory : public TypeBase
    {
    public:
        TypeInMemory(void);
        TypeInMemory(const TypeInMemory& other) = delete;
    };

    // Partial specialization of TypeInMemory for vectors.
    template<typename T>
    class TypeInMemory<std::vector<T>> : public TypeBase
    {
    public:
        TypeInMemory(void);
        TypeInMemory(const TypeInMemory& other) = delete;
    };


    class Group;
    class DataSetBase : public HDF5Object
    {
    protected:
        DataSetBase(std::string name)
          : HDF5Object(name)
        { }

        static std::string createName(const HDF5File::Group& group,
                                        std::string dsetName);

    public:
        DataSetBase(void) = delete;
        DataSetBase(const DataSetBase& other) = delete;

        ~DataSetBase(void)
        {
            H5Dclose(getId());
        }
    };

    // Templatized base class.
    // We do this so that we can implement the create and open 
    // constructors once rather than in each of the derived classes.
    template<typename T>
    class DataSetTBase : public DataSetBase
    {
    public:
        DataSetTBase(void) = delete;
        DataSetTBase(const DataSetTBase<T>& other) = delete;

        // Create data set.
        DataSetTBase(const HDF5File::Group& group,
                        std::string dsetName,
                        const DataSpace& dspace);

        // Open existing data set.
        DataSetTBase(const HDF5File::Group& group, std::string dsetName);
    };

    template<typename T>
    class DataSet : public DataSetTBase<T> {
    public:
        DataSet(void) = delete;
        DataSet(const DataSet& other) = delete;

        // Create data set.
        DataSet(const HDF5File::Group& group,
                std::string dsetName,
                const DataSpace& dspace)
          : DataSetTBase<T>(group, dspace, dsetName)
        { }

        // Open existing data set.
        DataSet(const HDF5File::Group& group, std::string dsetName)
          : DataSetTBase<T>(group, dsetName)
        { }

        void write(const T& data) const;
        T read(void) const;
    };

    // Partial specialization for vector of T.
    // At least as of C++11, we cannot do partial specialization of
    // a member function so we have to declare a whole class.
    template<typename T>
    class DataSet<std::vector<T>> : public DataSetTBase<std::vector<T>> {
    public:
        DataSet(void) = delete;
        DataSet(const DataSet& other) = delete;

        // Create data set.
        DataSet(const HDF5File::Group& group,
                std::string dsetName,
                const DataSpace& dspace)
          : DataSetTBase<std::vector<T>>(group, dsetName, dspace)
        { }

        // Open existing data set.
        DataSet(const HDF5File::Group& group, std::string dsetName)
          : DataSetTBase<std::vector<T>>(group, dsetName)
        { }

        void write(const std::vector<T>& data) const;
        std::vector<T> read(void) const;
    };



#if READY
    // Specialization of DataSet for writing vector of vectors.
    template<typename T>
    class VectorsDataSet : public DataSet<std::vector<T>>
    {
    public:
        VectorsDataSet(void) = delete;
        VectorsDataSet(const VectorsDataSet& other) = delete;

        // Create data set.
        VectorsDataSet(const HDF5File::Group& group,
                const DataSpace& dspace,
                std::string name)
          : DataSet<std::vector<T>>(group, dspace, name)
        { }

        // Open existing data set.
        VectorsDataSet(const HDF5File::Group& group, std::string name)
          : DataSet<std::vector<T>>(group, name)
        { }

        void write(const std::vector<std::vector<T>>& data) const;
        std::vector<std::vector<T>> read(void) const;
    };
#endif // READY

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


    class AttributeBase : public HDF5Object
    {
    private:
        const HDF5Object& target;

    protected:
        // Construct an AttributeBase without creating/opening an attribute.
        AttributeBase(const HDF5Object& _target, std::string attrName)
          : HDF5Object(attrName),
            target(_target)
        { }

        // Construct an AttributeBase by opening an existing attribute.
        AttributeBase(const HDF5Object& _target,
                        std::string attrName,
                        bool /* junk to change signature */);

    public:
        // Create an AttributeBase.
        AttributeBase(void) = delete;
        AttributeBase(const AttributeBase& other) = delete;

        ~AttributeBase(void)
        {
            if(id != H5I_INVALID_HID)
            {
                H5Aclose(id);
            }
        }

        void Delete(void)
        {
            H5Adelete(target.getId(), getName().c_str());
            id = H5I_INVALID_HID;
        }
    };

    template<typename T>
    class Attribute : public AttributeBase
    {
    public:
        Attribute(void) = delete;
        Attribute(const Attribute& other) = delete;

        // Create an Attribute of the given HDF5 object (group or file).
        Attribute(const HDF5Object& target,
                    std::string attrName,
                    const DataSpace& ds);

        // Open an Attribute of the given HDF5 object (group or file).
        Attribute(const HDF5Object& target, std::string attrName)
          : AttributeBase(target, attrName, true)
        { }

        // Set the attribute's value to the given value.
        void setTo(const T& value) const;

        // Get the attribute's current value.
        T get(void) const;
    };

    // Partial specialization of Attribute to support vectors.
    // TODO Is there a way to do this without duplicating so much of the
    // class declaration?
    template<typename T>
    class Attribute<std::vector<T>> : public AttributeBase
    {
    public:
        Attribute(void) = delete;
        Attribute(const Attribute& other) = delete;

        // Create an attribute on the given object.
        Attribute(const HDF5Object& target,
                    std::string attrName,
                    const DataSpace& ds);

        // Open an attribute of the given object.
        Attribute(const HDF5Object& target, std::string attrName)
          : AttributeBase(target, attrName, true)
        { }

        // Set the attribute's value to the given value.
        void setTo(const std::vector<T>& value) const;


        // Get the attribute's current value.
        std::vector<T> get(void) const;
    };

    // Partial specialization of Attribute to support vectors of vectors.
    // TODO Is there a way to do this without duplicating so much of the
    // class declaration?
    template<typename T>
    class Attribute<std::vector<std::vector<T>>> : public AttributeBase
    {
    public:
        Attribute(void) = delete;
        Attribute(const Attribute& other) = delete;

        // Create an attribute on the given object.
        Attribute(const HDF5Object& target,
                    std::string attrName,
                    const DataSpace& ds);

        // Open an attribute of the given object.
        Attribute(const HDF5Object& target, std::string attrName)
          : AttributeBase(target, attrName, true)
        { }

        // Set the attribute's value to the given value.
        void setTo(const std::vector<std::vector<T>>& value) const;

        // Get the attribute's current value.
        std::vector<std::vector<T>> get(void) const;
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

// Ensure we have implementations of template classes.
#include "xolotlCore/io/HDF5FileType.h"
#include "xolotlCore/io/HDF5FileAttribute.h"
#include "xolotlCore/io/HDF5FileDataSpace.h"
#include "xolotlCore/io/HDF5FileDataSet.h"

#endif // XCORE_HDF5FILE_H

