#include <sstream>
#include "xolotlCore/io/HDF5File.h"
#include "xolotlCore/io/HDF5Exception.h"

namespace xolotlCore
{

HDF5File::Group::Group(const HDF5Object& _location,
                        fs::path path,
                        bool create)
  : HDF5Object(path.string()),
    location(_location)
{
    if(create)
    {
        id = H5Gcreate(location.getId(),
                        name.c_str(),
                        H5P_DEFAULT,    // link creation property list
                        H5P_DEFAULT,    // group creation property list
                        H5P_DEFAULT);   // group access property list
    }
    else
    {
        id = H5Gopen(location.getId(),
                        name.c_str(),
                        H5P_DEFAULT);   // group access property list
    }
    if(id < 0)
    {
        std::ostringstream estr;
        estr << "Unable to create/open Group " << name
            << " at location " << location.getName();
        throw HDF5Exception(estr.str());
    }
}

} // namespace xolotlCore
