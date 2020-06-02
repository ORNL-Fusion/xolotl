#ifndef XCORE_HDF5FILE_ATTRIBUTE_H
#define XCORE_HDF5FILE_ATTRIBUTE_H

#include <iostream>
#include <sstream>
#include <boost/range/counting_range.hpp>
#include <xolotl/io/HDF5FileType.h>

namespace xolotl
{
namespace io
{

template<typename T>
HDF5File::Attribute<T>::Attribute(const HDF5Object& target,
		std::string attrName, const DataSpace& ds) :
		AttributeBase(target, attrName) {
	setId(
			H5Acreate(target.getId(), getName().c_str(),
					TypeInFile<T>().getId(), ds.getId(), H5P_DEFAULT,
					H5P_DEFAULT));
	if (getId() < 0) {
		std::ostringstream estr;
		estr << "Unable to create attribute " << getName() << " on group "
				<< target.getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
void HDF5File::Attribute<T>::setTo(const T& value) const {
	auto status = H5Awrite(getId(), TypeInMemory<T>().getId(), &value);
	if (status < 0) {
		std::ostringstream estr;
		estr << "Unable to set value of attribute " << getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
T HDF5File::Attribute<T>::get(void) const {
	T ret;
	auto status = H5Aread(getId(), TypeInMemory<T>().getId(), &ret);
	if (status < 0) {
		std::ostringstream estr;
		estr << "Unable to read value of attribute " << getName();
		throw HDF5Exception(estr.str());
	}
	return ret;
}

//----------------------------------------------------------------------------
// Partial specialization for vectors.

template<typename T>
HDF5File::Attribute<std::vector<T>>::Attribute(const HDF5Object& target,
		std::string attrName, const DataSpace& ds) :
		AttributeBase(target, attrName) {
	TypeInFile<T> ftype;
	setTo(
			H5Acreate(target.getId(), getName().c_str(), ftype.getId(),
					ds.getId(), H5P_DEFAULT, H5P_DEFAULT));
	if (getId() < 0) {
		std::ostringstream estr;
		estr << "Unable to create attribute " << getName() << " on group "
				<< target.getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
void HDF5File::Attribute<std::vector<T>>::setTo(
		const std::vector<T>& value) const {
	// Verify the dimensionality of the attribute.
	SimpleDataSpace<1> dspace(*this);
	auto nItems = dspace.getDims()[0];
	if (nItems != value.size()) {
		throw HDF5Exception(
				"Size mismatch when setting vector-valued attribute");
	}

	// Write the data to the attribute.
	TypeInMemory<T> memType;
	auto status = H5Awrite(getId(), memType.getId(), value.data());
	if (status < 0) {
		std::ostringstream estr;
		estr << "Failed to write variable length vector data to attribute "
				<< getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
std::vector<T> HDF5File::Attribute<std::vector<T>>::get(void) const {
	// Ensure we have space for the data.
	SimpleDataSpace<1> dspace(*this);
	auto nItems = dspace.getDims()[0];
	std::vector<T> ret(nItems);
	T* pdata = const_cast<T*>(ret.data());

	// Read the data from the file.
	TypeInMemory<T> memType;
	auto status = H5Aread(getId(), memType.getId(), pdata);
	if (status < 0) {
		std::ostringstream estr;
		estr << "Unable to read value of vector attribute " << getName();
		throw HDF5Exception(estr.str());
	}

	return ret;
}

//----------------------------------------------------------------------------
// Partial specialization for vector of vectors.

template<typename T>
HDF5File::Attribute<std::vector<std::vector<T>>>::Attribute(
		const HDF5Object& target,
		std::string attrName,
		const DataSpace& ds)
: AttributeBase(target, attrName)
{
	TypeInFile<std::vector<T>> ftype;
	setId(H5Acreate(target.getId(),
					getName().c_str(),
					ftype.getId(),
					ds.getId(),
					H5P_DEFAULT,
					H5P_DEFAULT));
	if(getId() < 0)
	{
		std::ostringstream estr;
		estr << "Unable to create attribute " << getName() << " on group " << target.getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
void HDF5File::Attribute<std::vector<std::vector<T>>>::setTo(
		const std::vector<std::vector<T>>& value) const
{
	// Define variable-length data we will write.
	auto nVectors = value.size();
	std::vector<hvl_t> vlData(nVectors);
	for(auto i : boost::counting_range<uint32_t>(0, nVectors))
	{
		vlData[i].p = const_cast<void*>(static_cast<const void*>(value[i].data()));
		vlData[i].len = value[i].size();
	}

	// Write the data to the attribute in the file.
	TypeInMemory<std::vector<T>> memType;
	auto status = H5Awrite(getId(),
			memType.getId(),
			vlData.data());
	if(status < 0)
	{
		std::ostringstream estr;
		estr << "Failed to write variable length vector data to attribute " << getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
std::vector<std::vector<T>> HDF5File::Attribute<std::vector<std::vector<T>>>::get(void) const
{
	// Determine the dataspace for our data set.
	SimpleDataSpace<1> dspace(*this);
	auto nVectors = dspace.getDims()[0];

	// Read the data from teh file.
	std::vector<hvl_t> vlData(nVectors);
	TypeInMemory<std::vector<T>> memType;
	auto status = H5Aread(getId(),
			memType.getId(),
			vlData.data());
	if(status < 0)
	{
		std::ostringstream estr;
		estr << "Failed to read variable length vector data from attribute " << getName();
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

#if READY
//----------------------------------------------------------------------------
// Specialization for vector of nDims-array of 2-array.

template<typename T>
HDF5File::ExtentsAttribute<T>::ExtentsAttribute( const HDF5Object& target,
		std::string attrName,
		const DataSpace& ds)
: AttributeBase(target, attrName)
{
	TypeInFile<T> ftype;
	id = H5Acreate(target.getId(),
			name.c_str(),
			ftype.getId(),
			ds.getId(),
			H5P_DEFAULT,
			H5P_DEFAULT);
	if(id < 0)
	{
		std::ostringstream estr;
		estr << "Unable to create attribute " << name << " on group " << target.getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
void
HDF5File::ExtentsAttribute<T>::setTo(
		const std::vector<std::vector<std::array<T, 2>>>& value) const
{
	// TODO do I convert it to a flat array?
	std::vector<T> flatData;
	for(auto i : boost::counting_range<uint32_t>(0, value.size()))
	{
		assert(value[i].size() == value[0].size());
		for(auto d : boost::counting_range<uint32_t>(0, value[i].size()))
		{
			std::copy(value[i][d].begin(), value[i][d].end(),
					std::back_inserter(flatData));
		}
	}

	// Write the data to the attribute in the file.
	TypeInMemory<T> memType;
	auto status = H5Awrite(id,
			memType.getId(),
			flatData.data());
	if(status < 0)
	{
		std::ostringstream estr;
		estr << "Failed to write data to attribute " << getName();
		throw HDF5Exception(estr.str());
	}
}

template<typename T>
std::vector<std::vector<std::array<T, 2>>>
HDF5File::ExtentsAttribute<T>::get(void) const
{
	// Determine the dataspace for our data set.
	SimpleDataSpace<3> dspace(*this);
	auto const& dims = dspace.getDims();

	// How many values total do we expect?
	auto nValues = std::accumulate(dims.begin(), dims.end(), 1,
			std::multiplies<hsize_t>());

	// Make space for the values.
	std::vector<T> flatData(nValues);

	// Read the data from the file.
	TypeInMemory<T> memType;
	auto status = H5Aread(id,
			memType.getId(),
			flatData.data());
	if(status < 0)
	{
		std::ostringstream estr;
		estr << "Failed to read data from attribute " << getName();
		throw HDF5Exception(estr.str());
	}

	// Convert to our output format.
	std::vector<std::vector<std::array<T, 2>>> ret(dims[0]);
	auto flatDataIter = flatData.begin();
	for(auto meshNum : boost::counting_range<uint32_t>(0, dims[0]))
	{
		ret[meshNum].resize(dims[1]);
		for(auto d : boost::counting_range<uint32_t>(0, dims[1]))
		{
			ret[meshNum][d][0] = *flatDataIter;
			++flatDataIter;
			ret[meshNum][d][1] = *flatDataIter;
			++flatDataIter;
		}
	}
	assert(flatDataIter == flatData.end());

	// Release memory allocated by the HDF5 library.
	H5Dvlen_reclaim(memType.getId(),
			dspace.getId(),
			H5P_DEFAULT,
			flatData.data());

	return ret;
}
#endif // READY

//----------------------------------------------------------------------------
// Full specialization for strings.

template<>
std::string
HDF5File::Attribute<std::string>::get(void) const;

template<>
void
HDF5File::Attribute<std::string>::setTo(const std::string& value) const;

//----------------------------------------------------------------------------
// Full specialization for vectors of strings.

template<>
HDF5File::Attribute<std::vector<std::string>>::Attribute(
		const HDF5Object& target, std::string attrName, const DataSpace& ds);

template<>
std::vector<std::string>
HDF5File::Attribute<std::vector<std::string>>::get(void) const;

template<>
void
HDF5File::Attribute<std::vector<std::string>>::setTo(
		const std::vector<std::string>& value) const;

} // namespace io
} // namespace xolotl

#endif // XCORE_HDF5FILE_ATTRIBUTE_H
