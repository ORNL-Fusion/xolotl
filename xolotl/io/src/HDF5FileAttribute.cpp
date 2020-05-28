#include <iostream>
#include <sstream>
#include <xolotl/io/HDF5File.h>

namespace xolotlCore {

HDF5File::AttributeBase::AttributeBase(const HDF5Object& _target,
		std::string _attrName, bool /* ignored */) :
		LocatedHDF5Object(_target, _attrName) {
	setId(H5Aopen(getLocation().getId(), getName().c_str(),
	H5P_DEFAULT));
	if (getId() < 0) {
		std::ostringstream estr;
		estr << "Unable to open attribute " << getName() << " on group "
				<< getLocation().getName();
		throw HDF5Exception(estr.str());
	}
}

//----------------------------------------------------------------------------
// Full specialization for string.

template<>
void HDF5File::Attribute<std::string>::setTo(const std::string& value) const {
	const char* pvalue = value.c_str();
	TypeInMemory<std::string> memType;
	auto status = H5Awrite(getId(), memType.getId(), &pvalue);
	if (status < 0) {
		std::ostringstream estr;
		estr << "Unable to set value of string attribute " << getName();
		throw HDF5Exception(estr.str());
	}
}

template<>
std::string HDF5File::Attribute<std::string>::get(void) const {
	std::string ret;

	char* pvalue = nullptr;
	TypeInMemory<std::string> memType;
	auto status = H5Aread(getId(), memType.getId(), &pvalue);
	if (status < 0) {
		std::ostringstream estr;
		estr << "Unable to get value of string attribute " << getName();
		throw HDF5Exception(estr.str());
	}
	ret = pvalue;
	free(pvalue);

	return ret;
}

//----------------------------------------------------------------------------
// Full specialization for vector of strings.

template<>
HDF5File::Attribute<std::vector<std::string> >::Attribute(
		const HDF5Object& target, std::string attrName, const DataSpace& ds) :
		AttributeBase(target, attrName) {
	TypeInFile<std::string> ftype;
	setId(
			H5Acreate(target.getId(), getName().c_str(), ftype.getId(),
					ds.getId(),
					H5P_DEFAULT,
					H5P_DEFAULT));
	if (getId() < 0) {
		std::ostringstream estr;
		estr << "Unable to create attribute " << getName() << " on group "
				<< target.getName();
		throw HDF5Exception(estr.str());
	}
}

template<>
std::vector<std::string> HDF5File::Attribute<std::vector<std::string> >::get(
		void) const {
	std::vector<std::string> ret;

	// Determine how many strings there are.
	SimpleDataSpace<1> dspace(*this);
	auto nItems = dspace.getDims()[0];

	if (nItems > 0) {
		// Make space for the return values.
		std::vector<char*> cdata(nItems, nullptr);

		TypeInMemory<std::string> memType;
		auto status = H5Aread(getId(), memType.getId(), cdata.data());
		if (status < 0) {
			std::ostringstream estr;
			estr << "Failed to read vector of strings from attribute "
					<< getName();
			throw HDF5Exception(estr.str());
		}

		// Convert to correct return format.
		ret.reserve(nItems);
		for (auto i : boost::counting_range<uint32_t>(0, nItems)) {
			// Save the current string.
			ret.emplace_back(cdata[i]);

			// Release memory allocated by HDF5 library.
			std::free(cdata[i]);
		}
	}

	return ret;
}

template<>
void HDF5File::Attribute<std::vector<std::string> >::setTo(
		const std::vector<std::string>& value) const {
	std::vector<const char*> cpdata(value.size());
	for (auto i : boost::counting_range<uint32_t>(0, value.size())) {
		cpdata[i] = value[i].c_str();
	}

	// Write the data to the attribute.
	TypeInMemory<std::string> memType;
	auto status = H5Awrite(getId(), memType.getId(), cpdata.data());
	if (status < 0) {
		std::ostringstream estr;
		estr << "Failed to write vector of strings to attribute " << getName();
		throw HDF5Exception(estr.str());
	}
}

} // namespace xolotlCore
