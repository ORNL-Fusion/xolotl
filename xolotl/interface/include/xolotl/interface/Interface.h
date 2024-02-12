#pragma once

#include <memory>

#include <mpi.h>

#include <xolotl/interface/IXolotlInterface.h>

namespace xolotl
{
namespace interface
{
std::unique_ptr<IXolotlInterface>
makeXolotlInterface(int& argc, const char* argv[]);
} /* namespace interface */
} /* namespace xolotl */
