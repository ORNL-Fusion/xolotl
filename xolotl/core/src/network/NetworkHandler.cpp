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
		std::cout << "\nLoaded network of size " << _network->getDOF() << "."
				  << std::endl;
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
