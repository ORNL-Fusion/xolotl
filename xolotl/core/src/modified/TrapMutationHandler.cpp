// Includes
#include <algorithm>
#include <iostream>

#include <mpi.h>

#include <xolotl/core/modified/TrapMutationHandler.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace modified
{
void
TrapMutationHandler::initialize(network::IReactionNetwork& network,
	network::IReactionNetwork::SparseFillMap& dfill, IdType nx, IdType ny,
	IdType nz)
{
	// This method fills two vectors to define the modified trap-mutation: for
	// the first one, the first value corresponds to the depth at which the He1
	// cluster undergo trap-mutation (if the value is negative it means that it
	// doesn't TM), the second value correspond to He2, etc.; the second vector
	// gives the size of the vacancies into which He trap-mutates. Information
	// about desorption is also initialized here.
	initializeDepthSize(network.getClusterCommon(0).getTemperature(0));

	// Add the needed reaction (dissociation) connectivity
	// Each (He_i)(V) cluster and I clusters are connected to He_i

	using NetworkType = network::IPSIReactionNetwork;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto specIdHe = psiNetwork->getHeliumSpeciesId();
	auto specIdI = psiNetwork->getInterstitialSpeciesId();
	auto specIdV = psiNetwork->getVacancySpeciesId();

	// Initialize the composition
	using AmountType = network::IReactionNetwork::AmountType;
	auto comp = std::vector<AmountType>(psiNetwork->getSpeciesListSize(), 0);

	// Loop on helium clusters from size 1 to 7
	for (std::size_t i = 1; i <= 7; i++) {
		// Skip if the depth is not defined
		if (depthVec[i - 1] < 0.0)
			continue;

		comp[specIdHe()] = i;
		comp[specIdI()] = 0;
		comp[specIdV()] = 0;
		auto heClusterId = psiNetwork->findClusterId(comp);
		// Check that the helium cluster is present in the network
		if (heClusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe helium cluster of size " +
				std::to_string(i) +
				"is not present in the network, "
				"cannot use the trap-mutation option!");
		}

		// The helium cluster is connected to itself
		// TODO: check if this create doublons and if this is a problem
		dfill[heClusterId].emplace_back(heClusterId);

		// Check the desorption
		if (i == desorp.size)
			desorp.id = heClusterId;

		// Get the size of the I/V for this helium size
		auto trapSize = sizeVec[i - 1];

		// Get the corresponding I cluster
		comp[specIdHe()] = 0;
		comp[specIdI()] = trapSize;
		auto iClusterId = psiNetwork->findClusterId(comp);
		// Check that the interstital cluster is present in the network
		if (iClusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe interstital cluster of size " +
				std::to_string(trapSize) +
				"is not present in the network, "
				"cannot use the trap-mutation option!");
		}

		// The interstitial cluster is connected to He
		if (!enableReducedJacobian)
			dfill[iClusterId].emplace_back(heClusterId);

		// Get the corresponding HeV cluster
		comp[specIdHe()] = i;
		comp[specIdV()] = trapSize;
		comp[specIdI()] = 0;
		auto heVClusterId = psiNetwork->findClusterId(comp);
		// Check that the HeV cluster is present in the network
		if (heVClusterId == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe HeV cluster of sizes " +
				std::to_string(i) + " " + std::to_string(trapSize) +
				"is not present in the network, "
				"cannot use the trap-mutation option!");
		}

		// The HeV cluster is connected to He
		if (!enableReducedJacobian)
			dfill[heVClusterId].emplace_back(heClusterId);
	}

	// Update the bubble bursting rate
	updateTrapMutationRate(network.getLargestRate());

	return;
}

void
TrapMutationHandler::initializeIndex1D(IdType surfacePos,
	network::IReactionNetwork& network,
	std::vector<advection::IAdvectionHandler*> advectionHandlers,
	std::vector<double> grid, IdType nx, IdType xs)
{
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	using NetworkType = network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;

	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto numSpecies = psiNetwork->getSpeciesListSize();
	auto specIdHe = psiNetwork->getHeliumSpeciesId();
	auto specIdI = psiNetwork->getInterstitialSpeciesId();
	auto specIdV = psiNetwork->getVacancySpeciesId();

	// No GB trap mutation handler in 1D for now

	// Create the temporary 2D vector
	auto& temp2DVector = tmBubbles.emplace_back();
	// Create the temporary 1D vector
	auto& temp1DVector = temp2DVector.emplace_back();

	// Loop on the grid points in the depth direction
	for (auto i = 0; i < nx; i++) {
		auto& indices = temp1DVector.emplace_back();
		// If we are on the left side of the surface there is no
		// modified trap-mutation
		if (i + xs <= surfacePos) {
			continue;
		}

		// Get the depth
		double depth =
			(grid[i + xs] + grid[i + xs + 1]) / 2.0 - grid[surfacePos + 1];
		double previousDepth =
			(grid[i + xs - 1] + grid[i + xs]) / 2.0 - grid[surfacePos + 1];

		// Loop on the depth vector
		for (auto l = 0; l < depthVec.size(); l++) {
			// Skip if the depth is negative
			if (depthVec[l] < 0.0)
				continue;
			// Check if a helium cluster undergo TM at this depth
			if (std::fabs(depth - depthVec[l]) < 0.01 ||
				(depthVec[l] - 0.01 < depth &&
					depthVec[l] - 0.01 > previousDepth)) {
				// Add the bubble of size l+1 to the indices
				auto comp = std::vector<AmountType>(numSpecies, 0);
				comp[specIdHe()] = l + 1;
				comp[specIdV()] = sizeVec[l];
				std::size_t heVClusterId = psiNetwork->findClusterId(comp);
				comp[specIdV()] = 0;
				std::size_t heClusterId = psiNetwork->findClusterId(comp);
				comp[specIdI()] = sizeVec[l];
				comp[specIdHe()] = 0;
				std::size_t iClusterId = psiNetwork->findClusterId(comp);
				indices.push_back({heVClusterId, heClusterId, iClusterId});
			}
		}
	}
}

void
TrapMutationHandler::initializeIndex2D(std::vector<IdType> surfacePos,
	network::IReactionNetwork& network,
	std::vector<advection::IAdvectionHandler*> advectionHandlers,
	std::vector<double> grid, IdType nx, IdType xs, IdType ny, double hy,
	IdType ys)
{
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	using NetworkType = network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto numSpecies = psiNetwork->getSpeciesListSize();
	auto specIdHe = psiNetwork->getHeliumSpeciesId();
	auto specIdI = psiNetwork->getInterstitialSpeciesId();
	auto specIdV = psiNetwork->getVacancySpeciesId();

	// Create a Sigma 3 trap mutation handler because it is the
	// only one available right now
	auto sigma3Handler = std::make_unique<Sigma3TrapMutationHandler>();
	auto sigma3DistanceVec = sigma3Handler->getDistanceVector();
	auto sigma3SizeVec = sigma3Handler->getSizeVector();

	// Create the temporary 2D vector
	auto& temp2DVector = tmBubbles.emplace_back();

	// Loop on the grid points in the Y direction
	temp2DVector.reserve(ny);
	for (auto j = 0; j < ny; j++) {
		// Create the temporary 1D vector
		auto& temp1DVector = temp2DVector.emplace_back();

		// Loop on the grid points in the depth direction
		for (auto i = 0; i < nx; i++) {
			// Create the list (vector) of indices at this grid point
			auto& indices = temp1DVector.emplace_back();

			// If we are on the left side of the surface there is no
			// modified trap-mutation
			if (i + xs <= surfacePos[j + ys]) {
				continue;
			}

			// Get the depth
			double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0 -
				grid[surfacePos[j + ys] + 1];
			double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0 -
				grid[surfacePos[j + ys] + 1];

			// Loop on the depth vector
			for (auto l = 0; l < depthVec.size(); l++) {
				// Skip if the depth is negative
				if (depthVec[l] < 0.0)
					continue;
				// Check if a helium cluster undergo TM at this depth
				if (std::fabs(depth - depthVec[l]) < 0.01 ||
					(depthVec[l] - 0.01 < depth &&
						depthVec[l] - 0.01 > previousDepth)) {
					// Add the bubble of size l+1 to the indices
					auto comp = std::vector<AmountType>(numSpecies, 0);
					comp[specIdHe()] = l + 1;
					comp[specIdV()] = sizeVec[l];
					auto heVClusterId = psiNetwork->findClusterId(comp);
					comp[specIdV()] = 0;
					auto heClusterId = psiNetwork->findClusterId(comp);
					comp[specIdI()] = sizeVec[l];
					comp[specIdHe()] = 0;
					auto iClusterId = psiNetwork->findClusterId(comp);
					indices.push_back({heVClusterId, heClusterId, iClusterId});
				}
			}

			// Get the Y position
			double yPos = (double)(j + ys) * hy;
			// Loop on the GB advection handlers
			for (auto n = 1; n < advectionHandlers.size(); n++) {
				// Get the location of the GB
				double location = advectionHandlers[n]->getLocation();
				// Get the current distance from the GB
				double distance = fabs(yPos - location);
				// Loop on the sigma 3 distance vector
				for (auto l = 0; l < sigma3DistanceVec.size(); l++) {
					// Skip if the distance is negative
					if (sigma3DistanceVec[l] < 0.0)
						continue;
					// Check if a helium cluster undergo TM at this depth
					if (std::fabs(distance - sigma3DistanceVec[l]) < 0.01) {
						// Add the bubble of size l+1 to the indices
						auto comp = std::vector<AmountType>(numSpecies, 0);
						comp[specIdHe()] = l + 1;
						comp[specIdV()] = sigma3SizeVec[l];
						auto heVClusterId = psiNetwork->findClusterId(comp);
						comp[specIdV()] = 0;
						auto heClusterId = psiNetwork->findClusterId(comp);
						comp[specIdI()] = sigma3SizeVec[l];
						comp[specIdHe()] = 0;
						auto iClusterId = psiNetwork->findClusterId(comp);
						std::array<std::size_t, 3> tempArray = {
							heVClusterId, heClusterId, iClusterId};
						// Check if this bubble is already
						// associated with this grid point.
						auto iter = std::find(
							indices.begin(), indices.end(), tempArray);
						if (iter == indices.end()) {
							// Add this bubble to the indices
							indices.push_back(tempArray);
						}
					}
				}
			}
		}
	}
}

void
TrapMutationHandler::initializeIndex3D(
	std::vector<std::vector<IdType>> surfacePos,
	network::IReactionNetwork& network,
	std::vector<advection::IAdvectionHandler*> advectionHandlers,
	std::vector<double> grid, IdType nx, IdType xs, IdType ny, double hy,
	IdType ys, IdType nz, double hz, IdType zs)
{
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	using NetworkType = network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	auto numSpecies = psiNetwork->getSpeciesListSize();
	auto specIdHe = psiNetwork->getHeliumSpeciesId();
	auto specIdI = psiNetwork->getInterstitialSpeciesId();
	auto specIdV = psiNetwork->getVacancySpeciesId();

	// Create a Sigma 3 trap mutation handler because it is the
	// only one available right now
	auto sigma3Handler = std::make_unique<Sigma3TrapMutationHandler>();
	auto sigma3DistanceVec = sigma3Handler->getDistanceVector();
	auto sigma3SizeVec = sigma3Handler->getSizeVector();

	// Loop on the grid points in the Z direction
	tmBubbles.reserve(nz);
	for (auto k = 0; k < nz; k++) {
		// Create the temporary 2D vector
		auto& temp2DVector = tmBubbles.emplace_back();

		// Loop on the grid points in the Y direction
		temp2DVector.reserve(ny);
		for (auto j = 0; j < ny; j++) {
			// Create the temporary 1D vector
			auto& temp1DVector = temp2DVector.emplace_back();

			// Loop on the grid points in the depth direction
			for (auto i = 0; i < nx; i++) {
				// Create the list (vector) of indices at this grid point
				auto& indices = temp1DVector.emplace_back();

				// If we are on the left side of the surface there is no
				// modified trap-mutation
				if (i + xs <= surfacePos[j + ys][k + zs]) {
					continue;
				}

				// Get the depth
				double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0 -
					grid[surfacePos[j + ys][k + zs] + 1];
				double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0 -
					grid[surfacePos[j + ys][k + zs] + 1];

				// Loop on the depth vector
				for (auto l = 0; l < depthVec.size(); l++) {
					// Skip if the depth is negative
					if (depthVec[l] < 0.0)
						continue;
					// Check if a helium cluster undergo TM at this depth
					if (std::fabs(depth - depthVec[l]) < 0.01 ||
						(depthVec[l] - 0.01 < depth &&
							depthVec[l] - 0.01 > previousDepth)) {
						// Add the bubble of size l+1 to the indices
						auto comp = std::vector<AmountType>(numSpecies, 0);
						comp[specIdHe()] = l + 1;
						comp[specIdV()] = sizeVec[l];
						auto heVClusterId = psiNetwork->findClusterId(comp);
						comp[specIdV()] = 0;
						auto heClusterId = psiNetwork->findClusterId(comp);
						comp[specIdI()] = sizeVec[l];
						comp[specIdHe()] = 0;
						auto iClusterId = psiNetwork->findClusterId(comp);
						indices.push_back(
							{heVClusterId, heClusterId, iClusterId});
					}
				}

				// Get the Y position
				double yPos = (double)(j + ys) * hy;
				// Loop on the GB advection handlers
				for (auto n = 1; n < advectionHandlers.size(); n++) {
					// Get the location of the GB
					double location = advectionHandlers[n]->getLocation();
					// Get the current distance from the GB
					double distance = fabs(yPos - location);

					// Loop on the sigma 3 distance vector
					for (auto l = 0; l < sigma3DistanceVec.size(); l++) {
						// Skip if the distance is negative
						if (sigma3DistanceVec[l] < 0.0)
							continue;
						// Check if a helium cluster undergo TM at this depth
						if (std::fabs(distance - sigma3DistanceVec[l]) < 0.01) {
							// Add the bubble of size l+1 to the indices
							auto comp = std::vector<AmountType>(numSpecies, 0);
							comp[specIdHe()] = l + 1;
							comp[specIdV()] = sigma3SizeVec[l];
							auto heVClusterId = psiNetwork->findClusterId(comp);
							comp[specIdV()] = 0;
							auto heClusterId = psiNetwork->findClusterId(comp);
							comp[specIdI()] = sigma3SizeVec[l];
							comp[specIdHe()] = 0;
							auto iClusterId = psiNetwork->findClusterId(comp);
							std::array<std::size_t, 3> tempArray = {
								heVClusterId, heClusterId, iClusterId};
							// Check if this bubble is already
							// associated with this grid point.
							auto iter = std::find(
								indices.begin(), indices.end(), tempArray);
							if (iter == indices.end()) {
								// Add this bubble to the indices
								indices.push_back(tempArray);
							}
						}
					}
				}
			}
		}
	}
}

void
TrapMutationHandler::updateTrapMutationRate(const double rate)
{
	// Multiply the biggest rate in the network by 1000.0
	// so that trap-mutation overcomes any other reaction
	kMutation = 1000.0 * rate;
}

void
TrapMutationHandler::updateDisappearingRate(double conc)
{
	// Set the rate to have an exponential decrease
	if (enableAttenuation) {
		kDis = exp(-4.0 * conc);
	}
}

void
TrapMutationHandler::computeTrapMutation(network::IReactionNetwork& network,
	double* concOffset, double* updatedConcOffset, IdType xi, IdType yj,
	IdType zk)
{
	// Initialize the rate of the reaction
	double rate = 0.0;

	// Loop on the list
	for (auto ids : tmBubbles[zk][yj][xi]) {
		auto bubbleIndex = std::get<0>(ids);
		auto heIndex = std::get<1>(ids);
		auto iIndex = std::get<2>(ids);

		// Get the initial concentration of helium
		double oldConc = concOffset[heIndex];

		// Check the desorption
		if (heIndex == desorp.id) {
			auto dof = network.getDOF();
			using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
				Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);
			// Get the left side rate (combination + emission)
			double totalRate = network.getLeftSideRate(dConcs, heIndex, xi + 1);
			// Define the trap-mutation rate taking into account the desorption
			rate = kDis * totalRate * (1.0 - desorp.portion) / desorp.portion;
		}
		else {
			rate = kDis * kMutation;
		}

		// Update the concentrations (the helium cluster loses its
		// concentration)
		updatedConcOffset[heIndex] -= rate * oldConc;
		updatedConcOffset[bubbleIndex] += rate * oldConc;
		updatedConcOffset[iIndex] += rate * oldConc;
	}
}

size_t
TrapMutationHandler::computePartialsForTrapMutation(
	network::IReactionNetwork& network, double* concOffset, double* val,
	IdType* indices, IdType xi, IdType yj, IdType zk)
{
	// Initialize the rate of the reaction
	double rate = 0.0;

	// Consider all bubbles at this grid point.
	// TODO Relying on convention for indices in indices/vals arrays is
	// error prone - could be done with multiple parallel arrays.
	uint32_t i = 0;
	for (auto ids : tmBubbles[zk][yj][xi]) {
		auto bubbleIndex = std::get<0>(ids);
		auto heIndex = std::get<1>(ids);
		auto iIndex = std::get<2>(ids);

		// Check the desorption
		if (heIndex == desorp.id) {
			auto dof = network.getDOF();
			using HostUnmanaged = Kokkos::View<double*, Kokkos::HostSpace,
				Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);
			// Get the left side rate (combination + emission)
			double totalRate = network.getLeftSideRate(dConcs, heIndex, xi + 1);
			// Define the trap-mutation rate taking into account the desorption
			rate = kDis * totalRate * (1.0 - desorp.portion) / desorp.portion;
		}
		else {
			rate = kDis * kMutation;
		}

		// Set the helium cluster partial derivative
		auto baseIndex = i * 3;
		indices[baseIndex] = heIndex;
		val[baseIndex] = -rate;

		if (!enableReducedJacobian) {
			// Set the bubble cluster partial derivative
			indices[(baseIndex) + 1] = bubbleIndex;
			val[(baseIndex) + 1] = rate;

			// Set the interstitial cluster partial derivative
			indices[(baseIndex) + 2] = iIndex;
			val[(baseIndex) + 2] = rate;
		}
		else {
			// Dummy numbers because it is off-diagonal
			indices[(baseIndex) + 1] = heIndex;
			val[(baseIndex) + 1] = 0.0;
			indices[(baseIndex) + 2] = heIndex;
			val[(baseIndex) + 2] = 0.0;
		}

		// Advance to next indices/vals index.
		++i;
	}

	return (tmBubbles[zk][yj][xi]).size();
}

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
