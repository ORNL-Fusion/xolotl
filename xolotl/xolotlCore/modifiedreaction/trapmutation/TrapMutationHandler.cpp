// Includes
#include <TrapMutationHandler.h>
#include <MathUtils.h>
#include <iostream>
#include <algorithm>
#include <mpi.h>

namespace xolotlCore {

void TrapMutationHandler::initialize(const IReactionNetwork& network, int nx,
		int xs, int ny, double hy, int ys, int nz, double hz, int zs) {
	// Add the needed reaction (dissociation) connectivity
	// Each (He_i)(V) cluster and I clusters are connected to He_i

	// Get the single interstitial cluster
	auto singleInterstitial = (IReactant *) network.get(Species::I, 1);
	// Get the double interstitial cluster
	auto doubleInterstitial = (IReactant *) network.get(Species::I, 2);
	// Get the triple interstitial cluster
	auto tripleInterstitial = (IReactant *) network.get(Species::I, 3);

	// Check He cluster is a type in the network
	auto& knownType = network.getKnownReactantTypes();

	// If the I clusters are not in the network,
	// there is no trap-mutation
	if (!singleInterstitial || !doubleInterstitial || !tripleInterstitial
			|| knownType.find(ReactantType::He) == knownType.end()) {
		// Clear the vector of HeV indices created by He undergoing trap-mutation
		// at each grid point
		tmBubbles.clear();

		// Change the value of ny and nz in 1D and 2D so that the same loop
		// works in every case
		if (nz == 0)
			nz = 1;
		if (ny == 0)
			ny = 1;

		// Loop on the grid points in the Z direction
		// TODO even with the reserve ops, this might involve lots of
		// separate memory allocations.
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
				temp1DVector.reserve(nx);
				for (int i = 0; i < nx; i++) {
					// Indicate no bubbles at this grid point.
					temp1DVector.emplace_back();
				}

				// Give the 1D vector to the 2D vector
				temp2DVector.emplace_back(temp1DVector);
			}

			// Give the 2D vector to the final vector
			tmBubbles.emplace_back(temp2DVector);
		}
		// Inform the user
		int procId;
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		if (procId == 0)
			std::cout
					<< "The modified trap-mutation won't happen because "
							"the interstitial clusters are missing or He type is missing."
					<< std::endl;

		return;
	}

	// Loop on the He clusters
	for (const auto& heMapItem : network.getAll(ReactantType::He)) {

		// Get the cluster and its size
		auto& cluster = static_cast<IReactant&>(*(heMapItem.second));
		int heSize = cluster.getSize();

		// The helium cluster is connected to itself
		cluster.setDissociationConnectivity(cluster.getId());

		// The single, double and triple interstitial clusters are connected to He
		singleInterstitial->setDissociationConnectivity(cluster.getId());
		doubleInterstitial->setDissociationConnectivity(cluster.getId());
		tripleInterstitial->setDissociationConnectivity(cluster.getId());

		// Loop on the bubbles
		for (const auto& hevMapItem : network.getAll(ReactantType::PSIMixed)) {

			// Get the bubble and its composition
			auto& bubble = static_cast<IReactant&>(*(hevMapItem.second));
			auto& comp = bubble.getComposition();

			// We are only interested in bubbles with one, two, or three vacancies
			if (comp[toCompIdx(Species::V)] > 3)
				continue;

			// Connect with He if the number of helium in the bubble is the same
			if (comp[toCompIdx(Species::He)] == heSize
					&& comp[toCompIdx(Species::D)] == 0
					&& comp[toCompIdx(Species::T)] == 0) {
				bubble.setDissociationConnectivity(cluster.getId());
			}
		}
	}

	// This method fills two vectors to define the modified trap-mutation: for the first one,
	// the first value corresponds to the depth at which the He1 cluster undergo trap-mutation
	// (if the value is negative it means that it doesn't TM), the second value correspond
	// to He2, etc.; the second vector gives the size of the vacancies into which He
	// trap-mutates. Information about desorption is also initialized here.
	initializeDepthSize(network.getTemperature());

	// Update the bubble bursting rate
	updateTrapMutationRate(network);

	return;
}

void TrapMutationHandler::initializeIndex1D(int surfacePos,
		const IReactionNetwork& network,
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int nx, int xs) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

	// No GB trap mutation handler in 1D for now

	// Create the temporary 2D vector
	ReactantRefVector2D temp2DVector;
	// Create the temporary 1D vector
	ReactantRefVector1D temp1DVector;

	// Loop on the grid points in the depth direction
	for (int i = 0; i < nx; i++) {
		// If we are on the left side of the surface there is no
		// modified trap-mutation
		if (i <= surfacePos) {
			temp1DVector.emplace_back();
			continue;
		}

		// Get the depth
		double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0
				- grid[surfacePos + 1];
		double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0
				- grid[surfacePos + 1];

		// Loop on the depth vector
		std::vector<std::reference_wrapper<IReactant> > indices;
		for (int l = 0; l < depthVec.size(); l++) {
			// Check if a helium cluster undergo TM at this depth
			if (std::fabs(depth - depthVec[l]) < 0.01
					|| (depthVec[l] - 0.01 < depth
							&& depthVec[l] - 0.01 > previousDepth)) {
				// Add the bubble of size l+1 to the indices
				// Loop on the bubbles
				for (auto const& heVMapItem : network.getAll(
						ReactantType::PSIMixed)) {
					// Get the bubble and its composition
					auto& bubble = static_cast<IReactant&>(*(heVMapItem.second));
					auto const& comp = bubble.getComposition();
					// Get the correct bubble
					if (comp[toCompIdx(Species::He)] == l + 1
							&& comp[toCompIdx(Species::V)] == sizeVec[l]
							&& comp[toCompIdx(Species::D)] == 0
							&& comp[toCompIdx(Species::T)] == 0) {
						// Add this bubble to the indices
						indices.emplace_back(bubble);
					}
				}
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
		const IReactionNetwork& network,
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

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
			std::vector<std::reference_wrapper<IReactant> > indices;

			// If we are on the left side of the surface there is no
			// modified trap-mutation
			if (i <= surfacePos[j]) {
				temp1DVector.push_back(indices);
				continue;
			}

			// Get the depth
			double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0
					- grid[surfacePos[j] + 1];
			double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0
					- grid[surfacePos[j] + 1];

			// Loop on the depth vector
			for (int l = 0; l < depthVec.size(); l++) {
				// Check if a helium cluster undergo TM at this depth
				if (std::fabs(depth - depthVec[l]) < 0.01
						|| (depthVec[l] - 0.01 < depth
								&& depthVec[l] - 0.01 > previousDepth)) {
					// Add the bubble of size l+1 to the indices
					// Loop on the bubbles
					for (auto const& heVMapItem : network.getAll(
							ReactantType::PSIMixed)) {
						// Get the bubble and its composition
						auto& bubble =
								static_cast<IReactant&>(*(heVMapItem.second));
						auto const& comp = bubble.getComposition();
						// Get the correct bubble
						if (comp[toCompIdx(Species::He)] == l + 1
								&& comp[toCompIdx(Species::V)] == sizeVec[l]
								&& comp[toCompIdx(Species::D)] == 0
								&& comp[toCompIdx(Species::T)] == 0) {
							// Add this bubble to the indices
							indices.emplace_back(bubble);
						}
					}
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
						// Loop on the bubbles
						for (auto const& heVMapItem : network.getAll(
								ReactantType::PSIMixed)) {
							// Get the bubble and its composition
							auto& bubble =
									static_cast<IReactant&>(*(heVMapItem.second));
							auto const& comp = bubble.getComposition();
							// Get the correct bubble
							if (comp[toCompIdx(Species::He)] == l + 1
									&& comp[toCompIdx(Species::V)]
											== sigma3SizeVec[l]
									&& comp[toCompIdx(Species::D)] == 0
									&& comp[toCompIdx(Species::T)] == 0) {
								// Check if this bubble is already 
								// associated with this grid point.
								auto biter =
										std::find_if(indices.begin(),
												indices.end(),
												[&bubble](const IReactant& testReactant) {
													return testReactant.getId() == bubble.getId();
												});
								if (biter == indices.end()) {
									// Add this bubble to the indices
									indices.emplace_back(bubble);
								}
							}
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
		const IReactionNetwork& network,
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int nx, int xs, int ny, double hy, int ys,
		int nz, double hz, int zs) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	tmBubbles.clear();

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
				std::vector<std::reference_wrapper<IReactant> > indices;

				// If we are on the left side of the surface there is no
				// modified trap-mutation
				if (i <= surfacePos[j][k]) {
					temp1DVector.emplace_back(indices);
					continue;
				}

				// Get the depth
				double depth = (grid[i + xs] + grid[i + xs + 1]) / 2.0
						- grid[surfacePos[j][k] + 1];
				double previousDepth = (grid[i + xs - 1] + grid[i + xs]) / 2.0
						- grid[surfacePos[j][k] + 1];

				// Loop on the depth vector
				for (int l = 0; l < depthVec.size(); l++) {
					// Check if a helium cluster undergo TM at this depth
					if (std::fabs(depth - depthVec[l]) < 0.01
							|| (depthVec[l] - 0.01 < depth
									&& depthVec[l] - 0.01 > previousDepth)) {
						// Add the bubble of size l+1 to the indices
						// Loop on the bubbles
						for (auto const& heVMapItem : network.getAll(
								ReactantType::PSIMixed)) {
							// Get the bubble and its composition
							auto& bubble =
									static_cast<IReactant&>(*(heVMapItem.second));
							auto const& comp = bubble.getComposition();
							// Get the correct bubble
							if (comp[toCompIdx(Species::He)] == l + 1
									&& comp[toCompIdx(Species::V)] == sizeVec[l]
									&& comp[toCompIdx(Species::D)] == 0
									&& comp[toCompIdx(Species::T)] == 0) {
								// Add this bubble to the indices
								indices.emplace_back(bubble);
							}
						}
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
							// Loop on the bubbles
							for (auto const& heVMapItem : network.getAll(
									ReactantType::PSIMixed)) {
								// Get the bubble and its composition
								auto& bubble =
										static_cast<IReactant&>(*(heVMapItem.second));
								auto const& comp = bubble.getComposition();
								// Get the correct bubble
								if (comp[toCompIdx(Species::He)] == l + 1
										&& comp[toCompIdx(Species::V)]
												== sigma3SizeVec[l]
										&& comp[toCompIdx(Species::D)] == 0
										&& comp[toCompIdx(Species::T)] == 0) {
									// Check if this bubble is already
									// associated with this grid point.
									auto biter =
											std::find_if(indices.begin(),
													indices.end(),
													[&bubble](const IReactant& testReactant) {
														return testReactant.getId() == bubble.getId();
													});
									if (biter == indices.end()) {
										// Add this bubble to the indices
										indices.emplace_back(bubble);
									}
								}
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

void TrapMutationHandler::updateTrapMutationRate(
		const IReactionNetwork& network) {
	// Multiply the biggest rate in the network by 1000.0
	// so that trap-mutation overcomes any other reaction
	kMutation = 1000.0 * network.getBiggestRate();

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

void TrapMutationHandler::computeTrapMutation(const IReactionNetwork& network,
		double *concOffset, double *updatedConcOffset, int xi, int yj, int zk) {

	// Initialize the rate of the reaction
	double rate = 0.0;

	// Loop on the list
	for (IReactant const& currReactant : tmBubbles[zk][yj][xi]) {

		// Get the stored bubble and its ID
		auto const& bubble = static_cast<IReactant const&>(currReactant);
		auto bubbleIndex = bubble.getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		// Note this composition has nonzero entries for both He and I,
		// so we can't use the network's get function that takes a composition.
		auto& comp = bubble.getComposition();
		auto heCluster = (IReactant *) network.get(Species::He,
				comp[toCompIdx(Species::He)]);
		auto heIndex = heCluster->getId() - 1;

		// Get the interstitial cluster with the same number of I as the number
		// of vacancies in the bubble and its ID
		auto iCluster = (IReactant *) network.get(Species::I,
				comp[toCompIdx(Species::V)]);
		auto iIndex = iCluster->getId() - 1;

		// Get the initial concentration of helium
		double oldConc = concOffset[heIndex];

		// Check the desorption
		if (comp[toCompIdx(Species::He)] == desorp.size) {
			// Get the left side rate (combination + emission)
			double totalRate = heCluster->getLeftSideRate(xi + 1);
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
		const IReactionNetwork& network, double *val, int *indices, int xi,
		int yj, int zk) {

	// Initialize the rate of the reaction
	double rate = 0.0;

	// Consider all bubbles at this grid point.
	// TODO Relying on convention for indices in indices/vals arrays is
	// error prone - could be done with multiple parallel arrays.
	uint32_t i = 0;
	for (IReactant const& currReactant : tmBubbles[zk][yj][xi]) {

		// Get the stored bubble and its ID
		auto const& bubble = static_cast<IReactant const&>(currReactant);
		auto bubbleIndex = bubble.getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		// Note this composition has non-zero entries for both He and I.
		// so we can't use the network's get function that takes a composition.
		auto const& comp = bubble.getComposition();
		auto heCluster = (IReactant *) network.get(Species::He,
				comp[toCompIdx(Species::He)]);
		auto heIndex = heCluster->getId() - 1;

		// Get the interstitial cluster with the same number of I as the number
		// of vacancies in the bubble and its ID
		auto iCluster = (IReactant *) network.get(Species::I,
				comp[toCompIdx(Species::V)]);
		auto iIndex = iCluster->getId() - 1;

		// Check the desorption
		if (comp[toCompIdx(Species::He)] == desorp.size) {
			// Get the left side rate (combination + emission)
			double totalRate = heCluster->getLeftSideRate(xi + 1);
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

