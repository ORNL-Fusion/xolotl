#include <fstream>
#include <iostream>
#include <stdexcept>

#include <xolotl/core/Constants.h>
#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
TemperatureHandler::~TemperatureHandler()
{
}

void
TemperatureHandler::initialize(const int dof)
{
	_dof = dof;
}
} // namespace temperature
} // namespace core
} // namespace xolotl
