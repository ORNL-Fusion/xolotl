// Includes
#include <TrapMutationHandler.h>
#include <experimental/PSIReactionNetwork.h>
#include <MathUtils.h>
#include <iostream>
#include <algorithm>
#include <mpi.h>

namespace xolotlCore {

void TrapMutationHandler::initialize(experimental::IReactionNetwork &network,
		xolotlCore::experimental::IReactionNetwork::SparseFillMap &dfill,
		int nx, int ny, int nz) {
	// This method fills two vectors to define the modified trap-mutation: for the first one,
	// the first value corresponds to the depth at which the He1 cluster undergo trap-mutation
	// (if the value is negative it means that it doesn't TM), the second value correspond
	// to He2, etc.; the second vector gives the size of the vacancies into which He
	// trap-mutates. Information about desorption is also initialized here.
	initializeDepthSize(network.getClusterCommon(0).getTemperature(0));

	// Add the needed reaction (dissociation) connectivity
	// Each (He_i)(V) cluster and I clusters are connected to He_i

	using NetworkType =
	experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// Initialize the composition
	NetworkType::Composition comp = NetworkType::Composition::zero();

	// Loop on helium clusters from size 1 to 7
	for (std::size_t i = 1; i <= 7; i++) {
		// Skip if the depth is not defined
		if (depthVec[i - 1] < 0.0)
			continue;

		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::He] = i;
		comp[NetworkType::Species::V] = 0;
		auto heCluster = psiNetwork->findCluster(comp, plsm::onHost);
		auto heClusterId = heCluster.getId();
		// Check that the helium cluster is present in the network
		if (heClusterId == NetworkType::invalidIndex()) {
			throw std::string(
					"\nThe helium cluster of size " + std::to_string(i)
							+ "is not present in the network, "
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
		comp[NetworkType::Species::He] = 0;
		comp[NetworkType::Species::I] = trapSize;
		auto iCluster = psiNetwork->findCluster(comp, plsm::onHost);
		auto iClusterId = iCluster.getId();
		// Check that the interstital cluster is present in the network
		if (iClusterId == NetworkType::invalidIndex()) {
			throw std::string(
					"\nThe interstital cluster of size "
							+ std::to_string(trapSize)
							+ "is not present in the network, "
									"cannot use the trap-mutation option!");
		}

		// The interstitial cluster is connected to He
		dfill[iClusterId].emplace_back(heClusterId);

		// Get the corresponding HeV cluster
		comp[NetworkType::Species::He] = i;
		comp[NetworkType::Species::V] = trapSize;
		comp[NetworkType::Species::I] = 0;
		auto heVCluster = psiNetwork->findCluster(comp, plsm::onHost);
		auto heVClusterId = heVCluster.getId();
		// Check that the HeV cluster is present in the network
		if (heVClusterId == NetworkType::invalidIndex()) {
			throw std::string(
					"\nThe HeV cluster of sizes " + std::to_string(i) + " "
							+ std::to_string(trapSize)
							+ "is not present in the network, "
									"cannot use the trap-mutation option!");
		}

		// The HeV cluster is connected to He
		dfill[heVClusterId].emplace_back(heClusterId);
	}

	// Update the bubble bursting rate
	updateTrapMutationRate(network.getLargestRate());

	return;
}

void TrapMutationHandler::initializeIndex1D(int surfacePos,
		experimental::IReactionNetwork &network,
		std::vector<IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	using NetworkType =
	experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// No GB trap mutation handler in 1D for now

	// Create the temporary 2D vector
	ReactantRefVector2D temp2DVector;
	// Create the temporary 1D vector
	ReactantRefVector1D temp1DVector;

	// Loop on the grid points in the depth direction
	for (int i = 0; i < nx; i++) {
		// If we are on the left side of the surface there is no
		// modified trap-mutation
		if (i + xs <= surfacePos) {
			temp1DVector.emplace_back();
			continue;
		}

		// Get the depth
		double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0
				- grid[surfacePos + 1];
		double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0
				- grid[surfacePos + 1];

		// Loop on the depth vector
		std::vector<std::tuple<std::size_t, std::size_t, std::size_t> > indices;
		for (int l = 0; l < depthVec.size(); l++) {
			// Check if a helium cluster undergo TM at this depth
			if (std::fabs(depth - depthVec[l]) < 0.01
					|| (depthVec[l] - 0.01 < depth
							&& depthVec[l] - 0.01 > previousDepth)) {
				// Add the bubble of size l+1 to the indices
				NetworkType::Composition comp;
				for (auto i : psiNetwork->getSpeciesRange()) {
					comp[i] = 0;
				}
				comp[NetworkType::Species::He] = l + 1;
				comp[NetworkType::Species::V] = sizeVec[l];
				auto heVCluster = psiNetwork->findCluster(comp, plsm::onHost);
				comp[NetworkType::Species::V] = 0;
				auto heCluster = psiNetwork->findCluster(comp, plsm::onHost);
				comp[NetworkType::Species::I] = sizeVec[l];
				comp[NetworkType::Species::He] = 0;
				auto iCluster = psiNetwork->findCluster(comp, plsm::onHost);
				indices.emplace_back(
						std::make_tuple(heVCluster.getId(), heCluster.getId(),
								iCluster.getId()));
			}
		}

		// Add indices to the index vector
		temp1DVector.emplace_back(indices);
	}

	// Give the 1D vector to the 2D vector
	temp2DVector.emplace_back(temp1DVector);

	// Give the 2D vector to the final vector
	tmBubbles.emplace_back(temp2DVector);

	return;
}

void TrapMutationHandler::initializeIndex2D(std::vector<int> surfacePos,
		experimental::IReactionNetwork &network,
		std::vector<IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	using NetworkType =
	experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// Create a Sigma 3 trap mutation handler because it is the
	// only one available right now
	auto sigma3Handler = new Sigma3TrapMutationHandler();
	auto sigma3DistanceVec = sigma3Handler->getDistanceVector();
	auto sigma3SizeVec = sigma3Handler->getSizeVector();

	// Create the temporary 2D vector
	ReactantRefVector2D temp2DVector;

	// Loop on the grid points in the Y direction
	temp2DVector.reserve(ny);
	for (int j = 0; j < ny; j++) {

		// Create the temporary 1D vector
		ReactantRefVector1D temp1DVector;

		// Loop on the grid points in the depth direction
		for (int i = 0; i < nx; i++) {
			// Create the list (vector) of indices at this grid point
			std::vector<std::tuple<std::size_t, std::size_t, std::size_t> > indices;

			// If we are on the left side of the surface there is no
			// modified trap-mutation
			if (i + xs <= surfacePos[j + ys]) {
				temp1DVector.push_back(indices);
				continue;
			}

			// Get the depth
			double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0
					- grid[surfacePos[j + ys] + 1];
			double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0
					- grid[surfacePos[j + ys] + 1];

			// Loop on the depth vector
			for (int l = 0; l < depthVec.size(); l++) {
				// Check if a helium cluster undergo TM at this depth
				if (std::fabs(depth - depthVec[l]) < 0.01
						|| (depthVec[l] - 0.01 < depth
								&& depthVec[l] - 0.01 > previousDepth)) {
					// Add the bubble of size l+1 to the indices
					NetworkType::Composition comp;
					for (auto i : psiNetwork->getSpeciesRange()) {
						comp[i] = 0;
					}
					comp[NetworkType::Species::He] = l + 1;
					comp[NetworkType::Species::V] = sizeVec[l];
					auto heVCluster = psiNetwork->findCluster(comp,
							plsm::onHost);
					comp[NetworkType::Species::V] = 0;
					auto heCluster = psiNetwork->findCluster(comp,
							plsm::onHost);
					comp[NetworkType::Species::I] = sizeVec[l];
					comp[NetworkType::Species::He] = 0;
					auto iCluster = psiNetwork->findCluster(comp, plsm::onHost);
					indices.emplace_back(
							std::make_tuple(heVCluster.getId(),
									heCluster.getId(), iCluster.getId()));
				}
			}

			// Get the Y position
			double yPos = (double) (j + ys) * hy;
			// Loop on the GB advection handlers
			for (int n = 1; n < advectionHandlers.size(); n++) {
				// Get the location of the GB
				double location = advectionHandlers[n]->getLocation();
				// Get the current distance from the GB
				double distance = fabs(yPos - location);
				// Loop on the sigma 3 distance vector
				for (int l = 0; l < sigma3DistanceVec.size(); l++) {
					// Check if a helium cluster undergo TM at this depth
					if (std::fabs(distance - sigma3DistanceVec[l]) < 0.01) {
						// Add the bubble of size l+1 to the indices
						NetworkType::Composition comp;
						for (auto i : psiNetwork->getSpeciesRange()) {
							comp[i] = 0;
						}
						comp[NetworkType::Species::He] = l + 1;
						comp[NetworkType::Species::V] = sigma3SizeVec[l];
						auto heVCluster = psiNetwork->findCluster(comp,
								plsm::onHost);
						comp[NetworkType::Species::V] = 0;
						auto heCluster = psiNetwork->findCluster(comp,
								plsm::onHost);
						comp[NetworkType::Species::I] = sizeVec[l];
						comp[NetworkType::Species::He] = 0;
						auto iCluster = psiNetwork->findCluster(comp,
								plsm::onHost);
						auto tempTuple = std::make_tuple(heVCluster.getId(),
								heCluster.getId(), iCluster.getId());
						// Check if this bubble is already
						// associated with this grid point.
						auto iter = std::find(indices.begin(), indices.end(),
								tempTuple);
						if (iter == indices.end()) {
							// Add this bubble to the indices
							indices.emplace_back(tempTuple);
						}
					}
				}
			}

			// Add indices to the index vector
			temp1DVector.push_back(indices);
		}

		// Give the 1D vector to the 2D vector
		temp2DVector.push_back(temp1DVector);
	}

	// Give the 2D vector to the final vector
	tmBubbles.push_back(temp2DVector);

	// Clear the memory
	delete sigma3Handler;

	return;
}

void TrapMutationHandler::initializeIndex3D(
		std::vector<std::vector<int> > surfacePos,
		experimental::IReactionNetwork &network,
		std::vector<IAdvectionHandler*> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys,
		int nz, double hz, int zs) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	using NetworkType =
	experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
	auto psiNetwork = dynamic_cast<NetworkType*>(&network);

	// Create a Sigma 3 trap mutation handler because it is the
	// only one available right now
	auto sigma3Handler = new Sigma3TrapMutationHandler();
	auto sigma3DistanceVec = sigma3Handler->getDistanceVector();
	auto sigma3SizeVec = sigma3Handler->getSizeVector();

	// Loop on the grid points in the Z direction
	tmBubbles.reserve(nz);
	for (int k = 0; k < nz; k++) {

		// Create the temporary 2D vector
		ReactantRefVector2D temp2DVector;

		// Loop on the grid points in the Y direction
		temp2DVector.reserve(ny);
		for (int j = 0; j < ny; j++) {
			// Create the temporary 1D vector
			ReactantRefVector1D temp1DVector;

			// Loop on the grid points in the depth direction
			for (int i = 0; i < nx; i++) {
				// Create the list (vector) of indices at this grid point
				std::vector<std::tuple<std::size_t, std::size_t, std::size_t> > indices;

				// If we are on the left side of the surface there is no
				// modified trap-mutation
				if (i + xs <= surfacePos[j + ys][k + zs]) {
					temp1DVector.emplace_back(indices);
					continue;
				}

				// Get the depth
				double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0
						- grid[surfacePos[j + ys][k + zs] + 1];
				double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0
						- grid[surfacePos[j + ys][k + zs] + 1];

				// Loop on the depth vector
				for (int l = 0; l < depthVec.size(); l++) {
					// Check if a helium cluster undergo TM at this depth
					if (std::fabs(depth - depthVec[l]) < 0.01
							|| (depthVec[l] - 0.01 < depth
									&& depthVec[l] - 0.01 > previousDepth)) {
						// Add the bubble of size l+1 to the indices
						NetworkType::Composition comp;
						for (auto i : psiNetwork->getSpeciesRange()) {
							comp[i] = 0;
						}
						comp[NetworkType::Species::He] = l + 1;
						comp[NetworkType::Species::V] = sizeVec[l];
						auto heVCluster = psiNetwork->findCluster(comp,
								plsm::onHost);
						comp[NetworkType::Species::V] = 0;
						auto heCluster = psiNetwork->findCluster(comp,
								plsm::onHost);
						comp[NetworkType::Species::I] = sizeVec[l];
						comp[NetworkType::Species::He] = 0;
						auto iCluster = psiNetwork->findCluster(comp,
								plsm::onHost);
						indices.emplace_back(
								std::make_tuple(heVCluster.getId(),
										heCluster.getId(), iCluster.getId()));
					}
				}

				// Get the Y position
				double yPos = (double) (j + ys) * hy;
				// Loop on the GB advection handlers
				for (int n = 1; n < advectionHandlers.size(); n++) {
					// Get the location of the GB
					double location = advectionHandlers[n]->getLocation();
					// Get the current distance from the GB
					double distance = fabs(yPos - location);

					// Loop on the sigma 3 distance vector
					for (int l = 0; l < sigma3DistanceVec.size(); l++) {
						// Check if a helium cluster undergo TM at this depth
						if (std::fabs(distance - sigma3DistanceVec[l]) < 0.01) {
							// Add the bubble of size l+1 to the indices
							NetworkType::Composition comp;
							for (auto i : psiNetwork->getSpeciesRange()) {
								comp[i] = 0;
							}
							comp[NetworkType::Species::He] = l + 1;
							comp[NetworkType::Species::V] = sigma3SizeVec[l];
							auto heVCluster = psiNetwork->findCluster(comp,
									plsm::onHost);
							comp[NetworkType::Species::V] = 0;
							auto heCluster = psiNetwork->findCluster(comp,
									plsm::onHost);
							comp[NetworkType::Species::I] = sizeVec[l];
							comp[NetworkType::Species::He] = 0;
							auto iCluster = psiNetwork->findCluster(comp,
									plsm::onHost);
							auto tempTuple = std::make_tuple(heVCluster.getId(),
									heCluster.getId(), iCluster.getId());
							// Check if this bubble is already
							// associated with this grid point.
							auto iter = std::find(indices.begin(),
									indices.end(), tempTuple);
							if (iter == indices.end()) {
								// Add this bubble to the indices
								indices.emplace_back(tempTuple);
							}
						}
					}
				}

				// Add indices to the index vector
				temp1DVector.push_back(indices);
			}

			// Give the 1D vector to the 2D vector
			temp2DVector.push_back(temp1DVector);
		}

		// Give the 2D vector to the final vector
		tmBubbles.push_back(temp2DVector);
	}

	// Clear the memory
	delete sigma3Handler;

	return;
}

void TrapMutationHandler::updateTrapMutationRate(const double rate) {
	// Multiply the biggest rate in the network by 1000.0
	// so that trap-mutation overcomes any other reaction
	kMutation = 1000.0 * rate;

	return;
}

void TrapMutationHandler::setAttenuation(bool isAttenuation) {
	attenuation = isAttenuation;

	return;
}

void TrapMutationHandler::updateDisappearingRate(double conc) {
	// Set the rate to have an exponential decrease
	if (attenuation)
		kDis = exp(-4.0 * conc);

	return;
}

void TrapMutationHandler::computeTrapMutation(
		experimental::IReactionNetwork &network, double *concOffset,
		double *updatedConcOffset, int xi, int yj, int zk) {

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
			using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);
			// Get the left side rate (combination + emission)
			double totalRate = network.getLeftSideRate(dConcs, heIndex, xi + 1);
			// Define the trap-mutation rate taking into account the desorption
			rate = kDis * totalRate * (1.0 - desorp.portion) / desorp.portion;
		} else {
			rate = kDis * kMutation;
		}

		// Update the concentrations (the helium cluster loses its concentration)
		updatedConcOffset[heIndex] -= rate * oldConc;
		updatedConcOffset[bubbleIndex] += rate * oldConc;
		updatedConcOffset[iIndex] += rate * oldConc;
	}

	return;
}

int TrapMutationHandler::computePartialsForTrapMutation(
		experimental::IReactionNetwork &network, double *concOffset,
		double *val, int *indices, int xi, int yj, int zk) {

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
			using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
			auto hConcs = HostUnmanaged(concOffset, dof);
			auto dConcs = Kokkos::View<double*>("Concentrations", dof);
			deep_copy(dConcs, hConcs);
			// Get the left side rate (combination + emission)
			double totalRate = network.getLeftSideRate(dConcs, heIndex, xi + 1);
			// Define the trap-mutation rate taking into account the desorption
			rate = kDis * totalRate * (1.0 - desorp.portion) / desorp.portion;
		} else {
			rate = kDis * kMutation;
		}

		// Set the helium cluster partial derivative
		auto baseIndex = i * 3;
		indices[baseIndex] = heIndex;
		val[baseIndex] = -rate;

		// Set the bubble cluster partial derivative
		indices[(baseIndex) + 1] = bubbleIndex;
		val[(baseIndex) + 1] = rate;

		// Set the interstitial cluster partial derivative
		indices[(baseIndex) + 2] = iIndex;
		val[(baseIndex) + 2] = rate;

		// Advance to next indices/vals index.
		++i;
	}

	return (tmBubbles[zk][yj][xi]).size();
}

}/* end namespace xolotlCore */

