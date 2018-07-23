#include <iostream>
#include <sstream>
#include "xolotlCore/io/HDF5File.h"
#include "xolotlCore/io/HDF5Exception.h"

namespace xolotlCore
{

HDF5File::ScalarDataSpace::ScalarDataSpace(void)
  : DataSpace(H5Screate(H5S_SCALAR))
{
    if(id < 0)
    {
        throw HDF5Exception("Failed to create scalar DataSpace");
    }
}

} // namespace xolotlCore

