#include <iostream>
#include <stdexcept>

#include <mpi.h>

#include <xolotl/core/network/NetworkHandler.h>

namespace xolotl
{
namespace core
{
namespace network
{
NetworkHandler::NetworkHandler(
	const options::IOptions& options, NetworkGeneratorFunction generatorFunc) :
	_network(generatorFunc(options))
{
	if (!_network) {
		throw std::runtime_error("Failed to load network");
	}

	int procId;
	MPI_Comm_rank(MPI_COMM_WORLD, &procId);
	if (procId == 0) {
		std::cout << "NetworkHandler: Loaded network of " << _network->getDOF()
				  << " DOF with: ";
		auto numSpecies = _network->getSpeciesListSize();
		for (auto id = SpeciesId(numSpecies); id; ++id) {
			auto speciesName = _network->getSpeciesName(id);
			std::cout << speciesName << " ";
		}
		std::cout << std::endl;
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
