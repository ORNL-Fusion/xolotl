#include <iostream>
#include "xolotlCore/io/HDF5File.h"
#include "xolotlCore/io/HDF5Exception.h"

namespace xolotlCore
{

template<>
HDF5File::TypeInFile<std::string>::TypeInFile(void)
  : HDF5File::TypeBase("std::string", H5Tcopy(H5T_C_S1), true)
{
    auto status = H5Tset_size(getId(), H5T_VARIABLE);
    if(status < 0)
    {
        throw HDF5Exception("Unable to construct variable-length string data type");
    }
}

template<>
HDF5File::TypeInMemory<std::string>::TypeInMemory(void)
  : HDF5File::TypeBase("std::string", H5Tcopy(H5T_C_S1), true)
{
    auto status = H5Tset_size(getId(), H5T_VARIABLE);
    if(status < 0)
    {
        throw HDF5Exception("Unable to construct variable-length string data type");
    }
}

} // namespace xolotlCore

