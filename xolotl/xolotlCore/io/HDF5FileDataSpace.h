#ifndef XCORE_HDF5FILE_DATASPACE_H
#define XCORE_HDF5FILE_DATASPACE_H

#include <iterator>
#include <sstream>
#include "xolotlCore/io/HDF5File.h"
#include "xolotlCore/io/HDF5Exception.h"

namespace xolotlCore
{

template<uint32_t Rank>
HDF5File::SimpleDataSpace<Rank>::SimpleDataSpace(const Dimensions& _dims)
  : DataSpace(H5Screate_simple(Rank, _dims.data(), nullptr)),
    dims(_dims)
{
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Failed to create simple data space for dims: ";
        std::copy(dims.begin(), dims.end(),
                std::ostream_iterator<uint32_t>(estr, " "));
        throw HDF5Exception(estr.str());
    }
}


template<uint32_t Rank>
HDF5File::SimpleDataSpace<Rank>::SimpleDataSpace(const Dimensions& _dims,
                                                const Dimensions& _maxDims)
  : DataSpace(H5Screate_simple(Rank, _dims.data(), _maxDims.data())),
    dims(_dims)
{
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Failed to create simple data space for dims: ";
        std::copy(dims.begin(), dims.end(),
                std::ostream_iterator<uint32_t>(estr, " "));
        throw HDF5Exception(estr.str());
    }
}


template<uint32_t Rank>
HDF5File::SimpleDataSpace<Rank>::SimpleDataSpace(const AttributeBase& attr)
  : DataSpace(H5Aget_space(attr.getId()))
{
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Failed to create simple data space for attribute " << attr.getName();
        throw HDF5Exception(estr.str());
    }

    // Retrieve the dimensions.
    auto nDims = H5Sget_simple_extent_dims(getId(),
                                        dims.data(),
                                        nullptr);
    if(nDims < 0)
    {
        std::ostringstream estr;
        estr << "Unable to obtain dimensions for data space for attribute " << attr.getName();
        throw HDF5Exception(estr.str());
    }
    assert(nDims == Rank);
}


template<uint32_t Rank>
HDF5File::SimpleDataSpace<Rank>::SimpleDataSpace(const DataSetBase& dsb)
  : DataSpace(H5Dget_space(dsb.getId()))
{
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Failed to construct simple data space for dataset " << dsb.getName();
        throw HDF5Exception(estr.str());
    }

    // Retrieve the dimensions.
    auto nDims = H5Sget_simple_extent_dims(getId(),
                                        dims.data(),
                                        nullptr);
    if(nDims < 0)
    {
        std::ostringstream estr;
        estr << "Unable to obtain dimensions for data space for dataset " << dsb.getName();
        throw HDF5Exception(estr.str());
    }
    assert(nDims == Rank);
}


template<uint32_t Rank>
void
HDF5File::SimpleDataSpace<Rank>::setDims(const Dimensions& _dims)
{
    auto status = H5Sset_extent_simple(getId(),
                                        Rank,
                                        _dims.data(),
                                        nullptr);
    if(status < 0)
    {
        std::ostringstream estr;
        estr << "Unable to update the dimensions for data space" << getName();
        throw HDF5Exception(estr.str());
    }
    dims = _dims;
}

} // namespace xolotlCore

#endif // XCORE_HDF5FILE_DATASPACE_H
