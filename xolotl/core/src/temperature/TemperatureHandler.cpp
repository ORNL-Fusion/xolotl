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
TemperatureHandler::initializeTemperature(const int dof,
	network::IReactionNetwork::SparseFillMap& ofillMap,
	network::IReactionNetwork::SparseFillMap& dfillMap,
	std::vector<double> grid)
{
	// Set dof
	_dof = dof;

	// Add the temperature to ofill
	ofillMap[_dof].emplace_back(_dof);

	// Add the temperature to dfill
	dfillMap[_dof].emplace_back(_dof);

	// keep the grd
	xGrid = grid;
}
} // namespace temperature
} // namespace core
} // namespace xolotl
