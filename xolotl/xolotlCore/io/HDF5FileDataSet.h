#ifndef XCORE_HDF5FILE_DATASET_H
#define XCORE_HDF5FILE_DATASET_H

#include <numeric>
#include "boost/range/counting_range.hpp"

namespace xolotlCore
{

inline
std::string
HDF5File::DataSetBase::createName(const HDF5File::Group& group,
                                std::string dsetName)
{
	std::ostringstream nstr;
	nstr << group.getName() << ':' << dsetName;
	return nstr.str();
}


template<typename T>
HDF5File::DataSet<T>::DataSet(const HDF5File::Group& group,
                                    const DataSpace& dspace,
                                    std::string dsetName)
  : DataSetBase(createName(group, dsetName))
{
    id = H5Dcreate(group.getId(),
                    dsetName.c_str(),
                    TypeInFile<T>().getId(),
                    dspace.getId(),
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT);
    if(id < 0)
    {
        std::ostringstream estr;
        estr << "Failed to create dataset " << getName();
        throw HDF5Exception(estr.str());
    }
}

template<typename T>
HDF5File::DataSet<T>::DataSet(const HDF5File::Group& group, std::string dsetName)
  : DataSetBase(createName(group, dsetName))
{
    id = H5Dopen(group.getId(),
                    dsetName.c_str(),
                    H5P_DEFAULT);
    if(id < 0)
    {
        std::ostringstream estr;
        estr << "Failed to open dataset " << getName();
        throw HDF5Exception(estr.str());
    }
}


template<typename T>
void
HDF5File::VectorsDataSet<T>::write(const std::vector<std::vector<T>>& data) const
{
    // Define variable-length data we will write.
    auto nDims = data.size();
    std::vector<hvl_t> vlData(nDims);
    for(auto i : boost::counting_range<uint32_t>(0, nDims))
    {
        vlData[i].p = const_cast<void*>(static_cast<const void*>(data[i].data()));
        vlData[i].len = data[i].size();
    }

    // Write the data to the file.
    TypeInMemory<std::vector<T>> memType;
    auto status = H5Dwrite(DataSet<std::vector<T>>::id,
                        memType.getId(),
                        H5S_ALL,
                        H5S_ALL,
                        H5P_DEFAULT,
                        vlData.data());
    if(status < 0)
    {
        throw HDF5Exception("Failed to write variable length vector data to dataset");
    }

    // We do *not* call H5vlen_reclaim here, 
    // because the data buffers were owned by the vectors passed
    // in by the caller.
}

template<typename T>
std::vector<std::vector<T>>
HDF5File::VectorsDataSet<T>::read(void) const
{
    // Determine the dataspace for our data set.
    // TODO is this really a dataspace with rank 1?
    SimpleDataSpace<1> dspace(*this);
    auto nVectors = dspace.getDims()[0];

    // Read the data from the file.
    std::vector<hvl_t> vlData(nVectors);
    TypeInMemory<std::vector<T>> memType;
    auto status = H5Dread(DataSet<std::vector<T>>::id,
                            memType.getId(),
                            H5S_ALL,
                            H5S_ALL,
                            H5P_DEFAULT,
                            vlData.data());
    if(status < 0)
    {
        std::ostringstream estr;
        estr << "Failed to read vector of vector data from dataset " 
            << this->getName();
        throw HDF5Exception(estr.str());
    }

    // Convert to our output format.
    std::vector<std::vector<T>> ret(nVectors);
    for(auto i : boost::counting_range<uint32_t>(0, nVectors))
    {
        for(auto j : boost::counting_range<uint32_t>(0, vlData[i].len))
        {
            T* currData = static_cast<T*>(vlData[i].p);
            ret[i].emplace_back(currData[j]);
        }
    }

    // Release memory allocated by the HDF5 library.
    H5Dvlen_reclaim(memType.getId(),
                    dspace.getId(),
                    H5P_DEFAULT,
                    vlData.data());

    return ret;
}

} // namespace xolotlCore

#endif // XCORE_HDF5FILE_DATASET_H
