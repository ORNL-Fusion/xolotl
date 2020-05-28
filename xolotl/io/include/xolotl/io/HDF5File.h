#ifndef XCORE_HDF5FILE_H
#define XCORE_HDF5FILE_H

#include <string>
#include <vector>
#include <array>
#include <mpi.h>
#include <xolotl/io/Filesystem.h>
#include <xolotl/io/HDF5Object.h>


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
          : HDF5Object("PropertyList", H5Pcreate(cls_id))
        { }

        /**
         * Release an open property list.
         */
        ~PropertyList(void)
        {
            H5Pclose(getId());
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
                H5Tclose(getId());
            }
        }
    };

    template<typename T>
    class TypeInFile : public TypeBase
    {
    private:
        hid_t BuildCompoundType(void) const;

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
    private:
        hid_t BuildCompoundType(void) const;

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


    class DataSetBase : public LocatedHDF5Object
    {
    protected:
        DataSetBase(const HDF5Object& loc, std::string name)
          : LocatedHDF5Object(loc, name)
        { }

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
        DataSetTBase(const HDF5Object& loc,
                        std::string dsetName,
                        const DataSpace& dspace);

        // Open existing data set.
        DataSetTBase(const HDF5Object& loc, std::string dsetName);
    };

    template<typename T>
    class DataSet : public DataSetTBase<T> {
    public:
        template<uint32_t dim0>
        using DataType2D = std::vector<std::array<T, dim0>>;

        DataSet(void) = delete;
        DataSet(const DataSet& other) = delete;

        // Create data set.
        DataSet(const HDF5Object& loc,
                std::string dsetName,
                const DataSpace& dspace)
          : DataSetTBase<T>(loc, dsetName, dspace)
        { }

        // Open existing data set.
        DataSet(const HDF5Object& loc, std::string dsetName)
          : DataSetTBase<T>(loc, dsetName)
        { }

        void write(const T& data) const;
        T read(void) const;

        template<uint32_t dim0>
        void parWrite2D(MPI_Comm comm,
                        uint32_t baseIdx,
                        const DataType2D<dim0>& data) const;
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
        DataSet(const HDF5Object& loc,
                std::string dsetName,
                const DataSpace& dspace)
          : DataSetTBase<std::vector<T>>(loc, dsetName, dspace)
        { }

        // Open existing data set.
        DataSet(const HDF5Object& loc, std::string dsetName)
          : DataSetTBase<std::vector<T>>(loc, dsetName)
        { }

        void write(const std::vector<T>& data) const;
        std::vector<T> read(void) const;
    };


    // Common base for all ragged data set classes.
    class RaggedDataSetBase {
    protected:    
        /// Concise name for type of flattened index metadata.
        using FlatStartingIndicesType = std::vector<uint32_t>;

        /// Suffix to add to data set names for starting indices dataset.
        static const std::string startIndicesDatasetNameSuffix;

        /// The MPI communicator used to access the file.
        MPI_Comm comm;


        /**
         * Create the data set.
         * Default and copy constructor explicitly disallowed.
         */
        RaggedDataSetBase(void) = delete;
        RaggedDataSetBase(const RaggedDataSetBase& other) = delete;

        /**
         * Create the data set.
         *
         * @param comm The MPI communicator used to access the file.
         */
        RaggedDataSetBase(MPI_Comm _comm)
          : comm(_comm)
        { }
    };


    // A DataSet for "Ragged" 2D data.  (I.e., 2D data where the 
    // number of items in dim 1 may vary).
    // Note that we do not inherit from DataSetTBase<std::vector<T>> because
    // we don't want to make our dataset have a variable-length type.
    template<typename T>
    class RaggedDataSet2D : public RaggedDataSetBase, public DataSetTBase<T> {
    public:
        /// Concise name for Ragged data type.
        using Ragged2DType = std::vector<std::vector<T>>;

    private:
        /// Concise name for type of flattened data.
        using FlatType = std::vector<T>;

        /**
         * Determine the number of values per grid point.
         *
         * @param data The ragged data set to be written.
         * @return Collection containing number of items per grid point.
         */
        static std::vector<uint32_t> findNumItemsByPoint(const Ragged2DType& data);

        /**
         * Build a dataspace for the flattened data representing the given
         * ragged data.
         *
         * @param _comm The MPI communicator used to access the file.
         * @param data The ragged data to be written.
         * @return A DataSpace describing the shape of the flattened dataset.
         */
        static std::unique_ptr<SimpleDataSpace<1>> buildDataSpace(
                                            MPI_Comm _comm,
                                            const Ragged2DType& data);

        /**
         * Read our part of the indexing metadata describing our
         * part of the flattened data set.
         *
         * @param baseX Index of the first X point we own.
         * @param numX Number of X points we own.
         * @return A vector containing the starting indices for 
         *              grid points we own
         */
        std::vector<uint32_t> readStartingIndices(int baseX, int numX) const;


        /**
         * Write our part of the indexing metadata describing our
         * part of the flattened data set.
         *
         * @param baseIdx Index of the first X point we own.
         * @param data Our part of the data to write.
         * @return Pair (globalBaseIdx, myNumItems) where 
         *              globalBaseIdx is index of our first item
         *              within the global flattened data set, and
         *              myNumItems is the number of items we own (and
         *              will write) from the flattened data set.
         */
        std::pair<uint32_t, uint32_t>
        writeStartingIndices(int baseX, const Ragged2DType& data) const;


        /**
         * Read our part of the ragged data set.
         *
         * @param globalStartingIndices Starting indices for each of
         * the grid points we own, plus one past so that we can compute
         * the total number of values we will read.
         * @return the part of the ragged data set that we own.
         */
        Ragged2DType readData(const std::vector<uint32_t>& globalStartingIndices) const ;

        /**
         * Write our part of the ragged data set.
         *
         * @param globalBaseIdx Index of first data item we own within
         *              the flattened data set.
         * @param myNumItems Number of items we own (and thus will write)
         *              within the flattened data set.
         * @param data The data to write.  (I.e., our part of the data. 
         */
        void writeData(uint32_t globalBaseIdx,
                        uint32_t myNumItems,
                        const Ragged2DType& data) const;

    public:
        RaggedDataSet2D(void) = delete;
        RaggedDataSet2D(const RaggedDataSet2D& other) = delete;

        /**
         * Create and write the data set.
         * Unlike the other DataSet classes that take a dataspace, 
         * this constructor takes the data itself so that it can
         * define the dataspaces for the flattened data and for
         * the indexing metadata that points into the flattened data.
         * Note that this constructor does *not* write the data, it
         * just defines the dataspaces for the data.
         *
         * @param comm The MPI communicator used to access the file.
         * @param loc The location (e.g., group) that contains our dataset.
         * @param dsetName The name of the dataset.
         * @param baseX Index of the first X point we own.
         * @param data The data to be written.
         */
        RaggedDataSet2D(MPI_Comm comm,
                        const HDF5Object& loc,
                        std::string dsetName,
                        int baseX,
                        const Ragged2DType& data);

        /**
         * Open an existing data set.
         *
         * @param comm The MPI communicator used to access the file.
         * @param loc The location (e.g., group) that contains our dataset.
         * @param dsetName The name of the dataset.
         */
        RaggedDataSet2D(MPI_Comm comm,
                        const HDF5Object& loc,
                        std::string dsetName);

        /**
         * Read data from an existing data set.
         *
         * @param baseX Index of the first X point we own.
         * @param numX Number of X points we own.
         * @return The data associated with X points in [baseX,baseX+numXs).
         */
        Ragged2DType read(int baseX, int numX) const;
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
        std::vector<std::vector<T>> read() const;
    };
#endif // READY

    // A group in the HDF5 file.
    class Group : public LocatedHDF5Object {
    public:                
        // Concise name for "ragged" 2D data.  (I.e., data where
        // number of items in dim 1 can vary.)
        template<typename T>
        using Ragged2DType = std::vector<std::vector<T>>;

        // Concise name for flattened starting index type used when
        // writing ragged n-D data.
        using FlatStartingIndicesType = std::vector<uint32_t>;

    private:
        /**
         * Write the starting indices for a 2D ragged dataset.
         *
         * @param comm The MPI communicator of our file.
         * @param commRank Our rank within the MPI communicator.
         * @param commSize The number of ranks in the MPI communicator.
         * @param baseX My rank's offset into dim 0 of the overall dataset.
         * @param data My part of the data to be written.
         */
        template<typename T>
        void writeStartingIndexDataset(MPI_Comm comm,
                                        int commRank,
                                        int commSize,
                                        int baseX,
                                        const Ragged2DType<T>& data) const;

        /**
         * Determine the starting indices of my items within
         * the flattened overall data set.
         *
         * @param comm The MPI communicator of our file.
         * @param commRank Our rank within the MPI communicator.
         * @param commSize The number of ranks in the MPI communicator.
         * @param The number of items associated with each of the points 
         *          that I own.
         */
        FlatStartingIndicesType
        findGlobalStartingIndices(MPI_Comm comm,
                            int commRank,
                            int commSize,
                            const FlatStartingIndicesType& localNumItems) const;

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
            H5Gclose(getId());
        }
    };


    class AttributeBase : public LocatedHDF5Object
    {
    protected:
        // Construct an AttributeBase without creating/opening an attribute.
        AttributeBase(const HDF5Object& _target, std::string attrName)
          : LocatedHDF5Object(_target, attrName)
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
            if(getId() != H5I_INVALID_HID)
            {
                H5Aclose(getId());
            }
        }

        void Delete(void)
        {
            H5Adelete(getLocation().getId(), getName().c_str());
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


    /**
     * The MPI communicator we are using.
     */
    MPI_Comm comm;


protected:
    /**
     * Determine whether we have a group at the named path in our file.
     *
     * @param path The path to check within our file.
     * @return True iff we have a group at the requested path.
     */
    bool hasGroup(fs::path path) const;

    /**
     * Access the MPI communicator we are using.
     *
     * @return The MPI communicator we are using.
     */
    MPI_Comm getComm(void) const  { return comm; }

	/**
	 * Close our file.
	 */
	void Close(void) {
		if(getId() != H5I_INVALID_HID) {
			H5Fclose(getId());
			setId(H5I_INVALID_HID);		
		}
	}

	/**
	 * Do the actual create/open of an HDF5 file.
	 *
	 * @param _path Path of file to create or open.
	 * @param _mode Access mode for creating/opening the file.
	 * @param _comm Communicator to use for accessing the file.
	 * @param par Whether to access the file using parallel I/O.
	 */
	void Open(fs::path _path,
				AccessMode _mode,
				MPI_Comm _comm,
				bool par);

public:
	/**
	 * Create or open an HDF5 file.
	 *
	 * @param _path Path of file to create or open.
	 * @param _mode Access mode for creating/opening the file.
	 * @param par Whether to access the file with parallel I/O.
	 */
	HDF5File(fs::path _path,
				AccessMode _mode,
				MPI_Comm _comm = MPI_COMM_WORLD,
				bool par = true )
	  : HDF5Object("/"),
		comm(_comm) {

		Open(_path, _mode, _comm, par);
	}

	/**
	 * Close the file if open and destroy the in-memory object.
	 */
	virtual ~HDF5File(void) {
		Close();
	}
};

} /* namespace xolotlCore */

// Ensure we have implementations of template classes.
#include <xolotl/io/HDF5FileType.h>
#include <xolotl/io/HDF5FileAttribute.h>
#include <xolotl/io/HDF5FileDataSpace.h>
#include <xolotl/io/HDF5FileDataSet.h>
#include <xolotl/io/HDF5FileGroup.h>

#endif // XCORE_HDF5FILE_H

