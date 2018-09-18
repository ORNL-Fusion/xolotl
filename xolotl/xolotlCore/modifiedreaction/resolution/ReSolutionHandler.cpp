// Includes
#include <ReSolutionHandler.h>
#include <NESuperCluster.h>
#include <algorithm>

namespace xolotlCore {

void ReSolutionHandler::initialize(const IReactionNetwork& network) {
	// Get the single xenon which connects to every cluster
	auto singleXenon = network.get(Species::Xe, 1);

	if (!singleXenon) {
		// Inform the user
		std::cout << "The re-solution won't happen because "
				"the single xenon cluster is missing." << std::endl;

		return;
	}

	// Loop on the clusters
	auto allClusters = network.getAll();
	std::for_each(allClusters.begin(), allClusters.end(),
			[&network,this,&singleXenon](IReactant& cluster) {
				// Get its size and radius
				auto size = cluster.getSize();
				auto radius = cluster.getReactionRadius();
				// Compute the re-soluted size
				double fraction = (0.949 * exp(-0.0703 * radius)
						+ (8.2326) / (1.0 + 7.982 * pow(radius, 2.0))
						* exp(-0.0371 * pow(radius, 2.0))) * 1.0e-4;
				int resolutedSize = (int) ((fraction * (double) size) + 0.5);

				// If the size is less than 1, the reaction should not happen
				if (resolutedSize > 0 && resolutedSize < size) {
					// Initialize the coefficients
					Array<double, 8> coefs;
					coefs.Init(0.0);
					// Different if it is a super cluster
					if (cluster.getType() == ReactantType::NESuper) {
						// Loop the bounds of the super cluster
						auto& superCluster = static_cast<NESuperCluster&> (cluster);
						int sectionWidth = superCluster.getSectionWidth();
						double dispersion = superCluster.getDispersion();
						double numXe = superCluster.getAverage();
						IReactant* previousSmaller = nullptr;
						for (int i = numXe + 1.0 - (double) sectionWidth / 2.0; i <= numXe + (double) sectionWidth / 2.0; i++) {
							// Update the quantities needed for the coefs
							double factor = ((double) i - numXe) / dispersion,
									distance = superCluster.getDistance(i), smallerFactor = 0.0;

							// Get the smaller cluster
							auto smallerCluster = network.get(Species::Xe, i - resolutedSize);
							// If it exists
							if (!smallerCluster) {
								// Get the super cluster
								smallerCluster = network.getSuperFromComp(i - resolutedSize, 0, 0, 0);
								auto smallerSuper = static_cast<NESuperCluster*> (smallerCluster);
								// Update the corresponding coefs
								smallerFactor = ((double) (i - resolutedSize) - smallerSuper->getAverage())
										/ smallerSuper->getDispersion();
							}


							// Check if it was the same as before
							if (smallerCluster != previousSmaller && previousSmaller) {
								// Set the connectivities
								cluster.setDissociationConnectivity(cluster.getId());
								cluster.setDissociationConnectivity(cluster.getMomentId());
								singleXenon->setDissociationConnectivity(cluster.getId());
								singleXenon->setDissociationConnectivity(cluster.getMomentId());
								previousSmaller->setDissociationConnectivity(cluster.getId());
								previousSmaller->setDissociationConnectivity(cluster.getMomentId());
								// Finish computing the factors
								coefs[2] = coefs[2] / ((double) sectionWidth);
								coefs[3] = coefs[3] / ((double) sectionWidth);
								coefs[4] = coefs[4] / ((double) sectionWidth);
								coefs[5] = coefs[5] / ((double) sectionWidth);
								coefs[6] = coefs[6] / ((double) sectionWidth);
								coefs[7] = coefs[7] / ((double) sectionWidth);
								// Add the size to the vector
								sizeVec.emplace_back(&cluster, previousSmaller, resolutedSize, coefs);

								// Reinitialize
								coefs.Init(0.0);
							}

							// Update the previous smaller cluster and the coefs
							previousSmaller = smallerCluster;
							coefs[0] += 1.0;
							coefs[1] += distance;
							coefs[2] += 1.0;
							coefs[3] += factor;
							coefs[4] += distance;
							coefs[5] += distance * factor;
							coefs[6] += smallerFactor;
							coefs[7] += distance * smallerFactor;
						}
						if (previousSmaller) {
							// Set the connectivities
							cluster.setDissociationConnectivity(cluster.getId());
							cluster.setDissociationConnectivity(cluster.getMomentId());
							singleXenon->setDissociationConnectivity(cluster.getId());
							singleXenon->setDissociationConnectivity(cluster.getMomentId());
							previousSmaller->setDissociationConnectivity(cluster.getId());
							previousSmaller->setDissociationConnectivity(cluster.getMomentId());
							// Finish computing the factors
							coefs[2] = coefs[2] / ((double) sectionWidth);
							coefs[3] = coefs[3] / ((double) sectionWidth);
							coefs[4] = coefs[4] / ((double) sectionWidth);
							coefs[5] = coefs[5] / ((double) sectionWidth);
							coefs[6] = coefs[6] / ((double) sectionWidth);
							coefs[7] = coefs[7] / ((double) sectionWidth);
							// Add the size to the vector
							sizeVec.emplace_back(&cluster, previousSmaller, resolutedSize, coefs);
						}
					}
					else {
						// Get the smaller cluster
						auto smallerCluster = network.get(Species::Xe, size - resolutedSize);

						// Set the connectivities
						cluster.setDissociationConnectivity(cluster.getId());
						singleXenon->setDissociationConnectivity(cluster.getId());
						smallerCluster->setDissociationConnectivity(cluster.getId());
						// Add the size to the vector
						coefs[0] = 1.0, coefs[2] = 1.0;
						sizeVec.emplace_back(&cluster, smallerCluster, resolutedSize, coefs);
					}
				}

			});

	return;
}

void ReSolutionHandler::updateReSolutionRate(double rate) {
	// Scale the fission rate with the standard one
	resolutionRate = 1.0e8 * 4.0 * rate; // The 4.0 factor is to get from the Xe rate to the fission rate

	return;
}

void ReSolutionHandler::computeReSolution(const IReactionNetwork& network,
		double *concOffset, double *updatedConcOffset, int xi, int xs, int yj,
		int zk) {
	// Get the single xenon which connects to every cluster
	auto singleXenon = network.get(Species::Xe, 1);
	int xenonId = singleXenon->getId() - 1;

	// Loop on the re-soluting clusters
	for (const auto& currPair : sizeVec) {
		// Get the larger cluster
		auto cluster = currPair.larger;
		int id = cluster->getId() - 1;
		int momId = cluster->getMomentId() - 1;
		// Get the re-solution cluster
		auto resoCluster = currPair.smaller;
		int resoId = resoCluster->getId() - 1;
		int resoMomId = resoCluster->getMomentId() - 1;

		// Get the initial concentration of the larger xenon cluster
		double l0 = concOffset[id], l1 = concOffset[momId];
		// Update the concentrations
		updatedConcOffset[id] -= resolutionRate
				* (currPair.coefs[2] * l0 + currPair.coefs[4] * l1);
		updatedConcOffset[momId] -= resolutionRate
				* (currPair.coefs[3] * l0 + currPair.coefs[5] * l1);
		updatedConcOffset[resoId] += resolutionRate
				* (currPair.coefs[2] * l0 + currPair.coefs[4] * l1);
		updatedConcOffset[resoMomId] += resolutionRate
				* (currPair.coefs[6] * l0 + currPair.coefs[7] * l1);
		updatedConcOffset[xenonId] += (double) currPair.size * resolutionRate
				* (currPair.coefs[0] * l0 + currPair.coefs[1] * l1);
	}

	return;
}

int ReSolutionHandler::computePartialsForReSolution(
		const IReactionNetwork& network, double *val, int *indices, int xi,
		int xs, int yj, int zk) {
	// Get the single xenon which connects to every cluster
	auto singleXenon = network.get(Species::Xe, 1);
	int xenonId = singleXenon->getId() - 1;

	// Loop on the re-soluting clusters
	int i = 0;
	for (const auto& currPair : sizeVec) {
		// Get the larger cluster
		auto cluster = currPair.larger;
		int id = cluster->getId() - 1;
		int momId = cluster->getMomentId() - 1;
		// Get the re-solution cluster
		auto resoCluster = currPair.smaller;
		int resoId = resoCluster->getId() - 1;
		int resoMomId = resoCluster->getMomentId() - 1;

		// Set the partial derivatives
		auto baseIndex = i * 10;
		indices[baseIndex] = id;
		val[baseIndex] = -resolutionRate * currPair.coefs[2];
		indices[(baseIndex) + 1] = momId;
		val[(baseIndex) + 1] = -resolutionRate * currPair.coefs[4];
		indices[(baseIndex) + 2] = id;
		val[(baseIndex) + 2] = -resolutionRate * currPair.coefs[3];
		indices[(baseIndex) + 3] = momId;
		val[(baseIndex) + 3] = -resolutionRate * currPair.coefs[5]; // Large cluster
		indices[(baseIndex) + 4] = resoId;
		val[(baseIndex) + 4] = resolutionRate * currPair.coefs[2];
		indices[(baseIndex) + 5] = resoMomId;
		val[(baseIndex) + 5] = resolutionRate * currPair.coefs[4];
		indices[(baseIndex) + 6] = resoId;
		val[(baseIndex) + 6] = resolutionRate * currPair.coefs[6];
		indices[(baseIndex) + 7] = resoMomId;
		val[(baseIndex) + 7] = resolutionRate * currPair.coefs[7]; // Smaller cluster
		indices[(baseIndex) + 8] = xenonId;
		val[(baseIndex) + 8] = (double) currPair.size * resolutionRate
				* currPair.coefs[0];
		indices[(baseIndex) + 9] = xenonId;
		val[(baseIndex) + 9] = (double) currPair.size * resolutionRate
				* currPair.coefs[1]; // Xe_1

		// Increment i
		i++;
	}

	return sizeVec.size();
}

}/* end namespace xolotlCore */

