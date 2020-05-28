#include <sstream>
#include <xolotl/io/HDF5File.h>
#include <xolotl/io/HDF5Exception.h>

namespace xolotlCore
{

HDF5File::Group::Group(const HDF5Object& _location,
                        fs::path path,
                        bool create)
  : LocatedHDF5Object(_location, path.string())
{
    if(create)
    {
        setId(H5Gcreate(getLocation().getId(),
                        getName().c_str(),
                        H5P_DEFAULT,    // link creation property list
                        H5P_DEFAULT,    // group creation property list
                        H5P_DEFAULT));   // group access property list
    }
    else
    {
        setId(H5Gopen(getLocation().getId(),
                        getName().c_str(),
                        H5P_DEFAULT));   // group access property list
    }
    if(getId() < 0)
    {
        std::ostringstream estr;
        estr << "Unable to create/open Group " << getName()
            << " at location " << getLocation().getName();
        throw HDF5Exception(estr.str());
    }
}

} // namespace xolotlCore
