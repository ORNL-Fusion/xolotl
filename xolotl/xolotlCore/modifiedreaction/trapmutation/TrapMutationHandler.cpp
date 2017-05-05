// Includes
#include <TrapMutationHandler.h>
#include <MathUtils.h>
#include <iostream>
#include <algorithm>

namespace xolotlCore {

void TrapMutationHandler::initialize(IReactionNetwork *network,
		std::vector<double> grid, int ny, double hy, int nz, double hz) {
	// Add the needed reaction (dissociation) connectivity
	// Each (He_i)(V) cluster and I clusters are connected to He_i

	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);
	// Get the single interstitial cluster
	auto singleInterstitial = (PSICluster *) network->get(iType, 1);
	// Get the double interstitial cluster
	auto doubleInterstitial = (PSICluster *) network->get(iType, 2);
	// Get the triple interstitial cluster
	auto tripleInterstitial = (PSICluster *) network->get(iType, 3);

	// If the I clusters are not in the network,
	// there is no trap-mutation
	if (!singleInterstitial || !doubleInterstitial || !tripleInterstitial) {
		// Clear the vector of HeV indices created by He undergoing trap-mutation
		// at each grid point
		indexVector.clear();

		// Change the value of ny and nz in 1D and 2D so that the same loop
		// works in every case
		if (nz == 0)
			nz = 1;
		if (ny == 0)
			ny = 1;

		// Loop on the grid points in the Z direction
		for (int k = 0; k < nz; k++) {
			// Create the temporary 2D vector
			std::vector<std::vector<std::vector<int> > > temp2DVector;
			// Loop on the grid points in the Y direction
			for (int j = 0; j < ny; j++) {
				// Create the temporary 1D vector
				std::vector<std::vector<int> > temp1DVector;
				// Loop on the grid points in the depth direction
				for (int i = 0; i < grid.size(); i++) {
					// Create the list (vector) of indices at this grid point
					std::vector<int> indices;
					// And give it empty to the 1D vector
					temp1DVector.push_back(indices);
				}

				// Give the 1D vector to the 2D vector
				temp2DVector.push_back(temp1DVector);
			}

			// Give the 2D vector to the final vector
			indexVector.push_back(temp2DVector);
		}
		// Inform the user
		std::cout << "The modified trap-mutation won't happen because "
				"the interstitial clusters are missing." << std::endl;

		return;
	}

	// Loop on the He clusters
	for (int i = 0; i < heClusters.size(); i++) {
		// Get the cluster and its size
		auto cluster = (PSICluster *) heClusters[i];
		int heSize = cluster->getSize();

		// The helium cluster is connected to itself
		cluster->setDissociationConnectivity(cluster->getId());

		// The single and double interstitial clusters are connected to He
		singleInterstitial->setDissociationConnectivity(cluster->getId());
		doubleInterstitial->setDissociationConnectivity(cluster->getId());
		tripleInterstitial->setDissociationConnectivity(cluster->getId());

		// Loop on the bubbles
		for (int j = 0; j < bubbles.size(); j++) {
			// Get the bubble and its composition
			auto bubble = (PSICluster *) bubbles[j];
			auto comp = bubble->getComposition();

			// We are only interested in bubbles with one, two, or three vacancies
			if (comp[vType] > 3)
				continue;

			// Connect with He if the number of helium in the bubble is the same
			if (comp[heType] == heSize) {
				bubble->setDissociationConnectivity(cluster->getId());
			}
		}
	}

	// This method fills two vectors to define the modified trap-mutation: for the first one,
	// the first value corresponds to the depth at which the He1 cluster undergo trap-mutation
	// (if the value is negative it means that it doesn't TM), the second value correspond
	// to He2, etc.; the second vector gives the size of the vacancies into which He
	// trap-mutates. Information about desorption is also initialized here.
	initializeDepthSize(network->getTemperature());

	// Update the bubble bursting rate
	updateTrapMutationRate(network);

	return;
}

void TrapMutationHandler::initializeIndex1D(int surfacePos,
		IReactionNetwork *network,
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	indexVector.clear();

	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);
	// No GB trap mutation handler in 1D for now

	// Create the temporary 2D vector
	std::vector<std::vector<std::vector<int> > > temp2DVector;
	// Create the temporary 1D vector
	std::vector<std::vector<int> > temp1DVector;

	// Loop on the grid points in the depth direction
	for (int i = 0; i < grid.size(); i++) {
		// Create the list (vector) of indices at this grid point
		std::vector<int> indices;

		// If we are on the left side of the surface there is no
		// modified trap-mutation
		if (i <= surfacePos) {
			temp1DVector.push_back(indices);
			continue;
		}

		// Get the depth
		double depth = grid[i] - grid[surfacePos];
		double previousDepth = grid[i - 1] - grid[surfacePos];

		// Loop on the depth vector
		for (int l = 0; l < depthVec.size(); l++) {
			// Check if a helium cluster undergo TM at this depth
			if (std::fabs(depth - depthVec[l]) < 0.01
					|| (depthVec[l] - 0.01 < depth
							&& depthVec[l] - 0.01 > previousDepth)) {
				// Add the bubble of size l+1 to the indices
				// Loop on the bubbles
				for (int m = 0; m < bubbles.size(); m++) {
					// Get the bubble and its composition
					auto bubble = (PSICluster *) bubbles[m];
					auto comp = bubble->getComposition();
					// Get the correct bubble
					if (comp[heType] == l + 1 && comp[vType] == sizeVec[l]) {
						// Add this bubble to the indices
						indices.push_back(m);
					}
				}
			}
		}

		// Add indices to the index vector
		temp1DVector.push_back(indices);
	}

	// Give the 1D vector to the 2D vector
	temp2DVector.push_back(temp1DVector);

	// Give the 2D vector to the final vector
	indexVector.push_back(temp2DVector);

	return;
}

void TrapMutationHandler::initializeIndex2D(std::vector<int> surfacePos,
		IReactionNetwork *network,
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int ny, double hy) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	indexVector.clear();

	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);
	// Create a Sigma 3 trap mutation handler because it is the
	// only one available right now
	auto sigma3Handler = new Sigma3TrapMutationHandler();
	auto sigma3DistanceVec = sigma3Handler->getDistanceVector();
	auto sigma3SizeVec = sigma3Handler->getSizeVector();

	// Create the temporary 2D vector
	std::vector<std::vector<std::vector<int> > > temp2DVector;
	// Loop on the grid points in the Y direction
	for (int j = 0; j < ny; j++) {
		// Create the temporary 1D vector
		std::vector<std::vector<int> > temp1DVector;
		// Loop on the grid points in the depth direction
		for (int i = 0; i < grid.size(); i++) {
			// Create the list (vector) of indices at this grid point
			std::vector<int> indices;

			// If we are on the left side of the surface there is no
			// modified trap-mutation
			if (i <= surfacePos[j]) {
				temp1DVector.push_back(indices);
				continue;
			}

			// Get the depth
			double depth = grid[i] - grid[surfacePos[j]];
			double previousDepth = grid[i - 1] - grid[surfacePos[j]];

			// Loop on the depth vector
			for (int l = 0; l < depthVec.size(); l++) {
				// Check if a helium cluster undergo TM at this depth
				if (std::fabs(depth - depthVec[l]) < 0.01
						|| (depthVec[l] - 0.01 < depth
								&& depthVec[l] - 0.01 > previousDepth)) {
					// Add the bubble of size l+1 to the indices
					// Loop on the bubbles
					for (int m = 0; m < bubbles.size(); m++) {
						// Get the bubble and its composition
						auto bubble = (PSICluster *) bubbles[m];
						auto comp = bubble->getComposition();
						// Get the correct bubble
						if (comp[heType] == l + 1
								&& comp[vType] == sizeVec[l]) {
							// Add this bubble to the indices
							indices.push_back(m);
						}
					}
				}
			}

			// Get the Y position
			double yPos = (double) j * hy;
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
						for (int m = 0; m < bubbles.size(); m++) {
							// Get the bubble and its composition
							auto bubble = (PSICluster *) bubbles[m];
							auto comp = bubble->getComposition();
							// Get the correct bubble
							if (comp[heType] == l + 1
									&& comp[vType] == sigma3SizeVec[l]) {
								// Check if this bubble is already in the indices
								if (std::find(indices.begin(), indices.end(), m)
										== indices.end()) {
									// Add this bubble to the indices
									indices.push_back(m);
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
	indexVector.push_back(temp2DVector);

	// Clear the memory
	delete sigma3Handler;

	return;
}

void TrapMutationHandler::initializeIndex3D(
		std::vector<std::vector<int> > surfacePos, IReactionNetwork *network,
		std::vector<IAdvectionHandler *> advectionHandlers,
		std::vector<double> grid, int ny, double hy, int nz, double hz) {
	// Clear the vector of HeV indices created by He undergoing trap-mutation
	// at each grid point
	indexVector.clear();

	// Get all the He clusters from the network
	auto heClusters = network->getAll(heType);
	// Get all the HeV bubbles from the network
	auto bubbles = network->getAll(heVType);
	// Create a Sigma 3 trap mutation handler because it is the
	// only one available right now
	auto sigma3Handler = new Sigma3TrapMutationHandler();
	auto sigma3DistanceVec = sigma3Handler->getDistanceVector();
	auto sigma3SizeVec = sigma3Handler->getSizeVector();

	// Loop on the grid points in the Z direction
	for (int k = 0; k < nz; k++) {
		// Create the temporary 2D vector
		std::vector<std::vector<std::vector<int> > > temp2DVector;
		// Loop on the grid points in the Y direction
		for (int j = 0; j < ny; j++) {
			// Create the temporary 1D vector
			std::vector<std::vector<int> > temp1DVector;
			// Loop on the grid points in the depth direction
			for (int i = 0; i < grid.size(); i++) {
				// Create the list (vector) of indices at this grid point
				std::vector<int> indices;

				// If we are on the left side of the surface there is no
				// modified trap-mutation
				if (i <= surfacePos[j][k]) {
					temp1DVector.push_back(indices);
					continue;
				}

				// Get the depth
				double depth = grid[i] - grid[surfacePos[j][k]];
				double previousDepth = grid[i - 1] - grid[surfacePos[j][k]];

				// Loop on the depth vector
				for (int l = 0; l < depthVec.size(); l++) {
					// Check if a helium cluster undergo TM at this depth
					if (std::fabs(depth - depthVec[l]) < 0.01
							|| (depthVec[l] - 0.01 < depth
									&& depthVec[l] - 0.01 > previousDepth)) {
						// Add the bubble of size l+1 to the indices
						// Loop on the bubbles
						for (int m = 0; m < bubbles.size(); m++) {
							// Get the bubble and its composition
							auto bubble = (PSICluster *) bubbles[m];
							auto comp = bubble->getComposition();
							// Get the correct bubble
							if (comp[heType] == l + 1
									&& comp[vType] == sizeVec[l]) {
								// Add this bubble to the indices
								indices.push_back(m);
							}
						}
					}
				}

				// Get the Y position
				double yPos = (double) j * hy;
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
							for (int m = 0; m < bubbles.size(); m++) {
								// Get the bubble and its composition
								auto bubble = (PSICluster *) bubbles[m];
								auto comp = bubble->getComposition();
								// Get the correct bubble
								if (comp[heType] == l + 1
										&& comp[vType] == sigma3SizeVec[l]) {
									// Check if this bubble is already in the indices
									if (std::find(indices.begin(),
											indices.end(), m)
											== indices.end()) {
										// Add this bubble to the indices
										indices.push_back(m);
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
		indexVector.push_back(temp2DVector);
	}

	// Clear the memory
	delete sigma3Handler;

	return;
}

void TrapMutationHandler::updateTrapMutationRate(IReactionNetwork *network) {
	// Multiply the biggest rate in the network by 1000.0
	// so that trap-mutation overcomes any other reaction
	kMutation = 1000.0 * network->getBiggestRate();

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

void TrapMutationHandler::computeTrapMutation(IReactionNetwork *network,
		double *concOffset, double *updatedConcOffset, int xi, int yj, int zk) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);
	// Initialyze the pointers to interstitial and helium clusters and their ID
	PSICluster * iCluster = nullptr, *heCluster = nullptr, *bubble = nullptr;
	int iIndex = -1, heIndex = -1, bubbleIndex = -1;

	// Initialize the rate of the reaction
	double rate = 0.0;

	// Get the pointer to list of indices at this grid point
	std::vector<int> * indices = &indexVector[zk][yj][xi];

	// Loop on the list
	for (int i = 0; i < indices->size(); i++) {
		// Get the stored bubble and its ID
		bubble = (PSICluster *) bubbles[indices->at(i)];
		bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		heCluster = (PSICluster *) network->get(heType, comp[heType]);
		heIndex = heCluster->getId() - 1;

		// Get the interstitial cluster with the same number of I as the number
		// of vacancies in the bubble and its ID
		iCluster = (PSICluster *) network->get(iType, comp[vType]);
		iIndex = iCluster->getId() - 1;

		// Get the initial concentration of helium
		double oldConc = concOffset[heIndex];

		// Check the desorption
		if (comp[heType] == desorp.size) {
			// Get the left side rate (combination + emission)
			double totalRate = heCluster->getLeftSideRate();
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
		IReactionNetwork *network, double *val, int *indices, int xi, int yj,
		int zk) {
	// Get all the HeV bubbles
	auto bubbles = network->getAll(heVType);
	// Initialyze the pointers to interstitial and helium clusters and their ID
	PSICluster * iCluster = nullptr, *heCluster = nullptr, *bubble = nullptr;
	int iIndex = -1, heIndex = -1, bubbleIndex = -1;

	// Initialize the rate of the reaction
	double rate = 0.0;

	// Get the pointer to list of indices at this grid point
	std::vector<int> * clusterIndices = &indexVector[zk][yj][xi];
	// Loop on the list
	for (int i = 0; i < clusterIndices->size(); i++) {
		// Get the stored bubble and its ID
		bubble = (PSICluster *) bubbles[clusterIndices->at(i)];
		bubbleIndex = bubble->getId() - 1;

		// Get the helium cluster with the same number of He and its ID
		auto comp = bubble->getComposition();
		heCluster = (PSICluster *) network->get(heType, comp[heType]);
		heIndex = heCluster->getId() - 1;

		// Get the interstitial cluster with the same number of I as the number
		// of vacancies in the bubble and its ID
		iCluster = (PSICluster *) network->get(iType, comp[vType]);
		iIndex = iCluster->getId() - 1;

		// Check the desorption
		if (comp[heType] == desorp.size) {
			// Get the left side rate (combination + emission)
			double totalRate = heCluster->getLeftSideRate();
			// Define the trap-mutation rate taking into account the desorption
			rate = kDis * totalRate * (1.0 - desorp.portion) / desorp.portion;
		} else {
			rate = kDis * kMutation;
		}

		// Set the helium cluster partial derivative
		indices[i * 3] = heIndex;
		val[i * 3] = -rate;

		// Set the bubble cluster partial derivative
		indices[(i * 3) + 1] = bubbleIndex;
		val[(i * 3) + 1] = rate;

		// Set the interstitial cluster partial derivative
		indices[(i * 3) + 2] = iIndex;
		val[(i * 3) + 2] = rate;
	}

	return clusterIndices->size();
}

}/* end namespace xolotlCore */

