// Includes
#include <ReSolutionHandler.h>
#include <experimental/NEReactionNetwork.h>
#include <algorithm>
#include <fstream>
#include <iostream>

namespace xolotlCore {

void ReSolutionHandler::initialize(experimental::IReactionNetwork &network,
		experimental::IReactionNetwork::SparseFillMap &dfill,
		double electronicStoppingPower) {

	using NetworkType =
	experimental::NEReactionNetwork;
	auto neNetwork = dynamic_cast<NetworkType*>(&network);
	// Get the single xenon which connects to every cluster
	NetworkType::Composition comp;
	for (auto i : neNetwork->getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[NetworkType::Species::Xe] = 1;
	auto singleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto xenonId = singleXenon.getId();

	if (xenonId == plsm::invalid<std::size_t>) {
		// Inform the user
		std::cout << "The re-solution won't happen because "
				"the single xenon cluster is missing." << std::endl;

		return;
	}

	// Set the fit variables depending on the electronic stopping power
	double y0 = 0.0, a1 = 0.0, b1 = 0.0, b2 = 0.0, c = 0.0;
	if (electronicStoppingPower > 0.87) {
		y0 = 11.0851, a1 = 1.5052, b1 = 0.0362, b2 = 0.0203, c = 3.4123;
	} else if (electronicStoppingPower > 0.82) {
		y0 = 10.6297, a1 = 1.3479, b1 = 0.0438, b2 = 0.0241, c = 4.2214;
	} else if (electronicStoppingPower > 0.77) {
		y0 = 10.1521, a1 = 1.1986, b1 = 0.0546, b2 = 0.0299, c = 5.4612;
	} else if (electronicStoppingPower > 0.71) {
		y0 = 9.1816, a1 = 0.949, b1 = 0.0703, b2 = 0.0371, c = 7.982;
	} else if (electronicStoppingPower > 0.67) {
		y0 = 8.6745, a1 = 0.8401, b1 = 0.0792, b2 = 0.0407, c = 9.6585;
	} else if (electronicStoppingPower > 0.62) {
		y0 = 7.6984, a1 = 0.6721, b1 = 0.1028, b2 = 0.0526, c = 14.272;
	} else if (electronicStoppingPower > 0.57) {
		y0 = 6.3925, a1 = 0.5025, b1 = 0.1411, b2 = 0.0727, c = 23.1967;
	} else if (electronicStoppingPower > 0.52) {
		y0 = 4.6175, a1 = 0.3433, b1 = 0.2284, b2 = 0.1276, c = 45.6624;
	} else {
		y0 = 2.3061, a1 = 0.2786, b1 = 1.1008, b2 = 1.605, c = 150.6689;
	}

	// Loop on the clusters
	constexpr auto speciesRange = NetworkType::getSpeciesRange();
	for (std::size_t id = 0; id < network.getNumClusters(); ++id) {
		auto cluster = neNetwork->getCluster(id, plsm::onHost);
		const auto &clReg = cluster.getRegion();
		NetworkType::Composition lo = clReg.getOrigin();
		NetworkType::Composition hi = clReg.getUpperLimitPoint();

		// Get its size
		auto size = lo[NetworkType::Species::Xe];

		// The re-soluted size is always 1
		int resolutedSize = 1;

		// If the resolutedSize is less than 1, or if the cluster is too small the reaction should not happen
		if (resolutedSize > 0 && resolutedSize < size && size >= minSize) {
			// Initialize the coefficients
			Array<double, 10> coefs;
			coefs.Init(0.0);
			// Different if it is a super cluster
			if (clReg.volume() > 1) {
				std::size_t sectionWidth = hi[NetworkType::Species::Xe]
						- lo[NetworkType::Species::Xe] - 1,
						smallerSectionWidth = 0.0;
				const auto &dispersion = clReg.dispersion();
				double numXe = (double) (hi[NetworkType::Species::Xe]
						+ lo[NetworkType::Species::Xe] - 1) / 2.0;
				auto previousSmallerId = plsm::invalid<std::size_t>;

				// Loop on all sizes
				for (int i = lo[NetworkType::Species::Xe];
						i < hi[NetworkType::Species::Xe]; i++) {
					// Update the quantities needed for the coefs
					double factor = ((double) i - numXe) / dispersion[0],
							smallerFactor = 0.0;

					// Get the smaller cluster
					comp[NetworkType::Species::Xe] = i - resolutedSize;
					// TODO: does this work for super clusters as well?
					auto smallerCluster = neNetwork->findCluster(comp,
							plsm::onHost);
					const auto &smallerReg = smallerCluster.getRegion();
					// If it exists and is super
					if (smallerReg.volume() > 1) {
						NetworkType::Composition smallerLo =
								smallerReg.getOrigin();
						NetworkType::Composition smallerHi =
								smallerReg.getUpperLimitPoint();
						// Update the corresponding coefs
						smallerFactor = ((double) (i - resolutedSize)
								- (double) (smallerHi[NetworkType::Species::Xe]
										+ smallerLo[NetworkType::Species::Xe]
										- 1) / 2.0)
								/ smallerReg.dispersion()[0];
					}

					// Check if it was the same as before
					if (smallerCluster.getId() != previousSmallerId
							&& previousSmallerId
									!= plsm::invalid<std::size_t>) {
						// Set the connectivities
						dfill[id].emplace_back(id);
						if (clReg.volume() > 1) {
							const auto &momIds = cluster.getMomentIds();
							for (auto j : speciesRange) {
								dfill[id].emplace_back(momIds[j()]);
								dfill[momIds[j()]].emplace_back(id);
								for (auto k : speciesRange) {
									dfill[momIds[j()]].emplace_back(
											momIds[k()]);
								}
							}
						}
						dfill[xenonId].emplace_back(id);
						if (clReg.volume() > 1) {
							const auto &momIds = cluster.getMomentIds();
							for (auto j : speciesRange) {
								dfill[xenonId].emplace_back(momIds[j()]);
							}
						}
						const auto &previousSmallerReg = neNetwork->getCluster(
								previousSmallerId, plsm::onHost).getRegion();
						dfill[previousSmallerId].emplace_back(id);
						if (clReg.volume() > 1) {
							const auto &momIds = cluster.getMomentIds();
							for (auto j : speciesRange) {
								dfill[previousSmallerId].emplace_back(
										momIds[j()]);
								if (previousSmallerReg.volume() > 1) {
									const auto &smallerMomIds =
											neNetwork->getCluster(
													previousSmallerId,
													plsm::onHost).getMomentIds();
									dfill[smallerMomIds[j()]].emplace_back(id);
									for (auto k : speciesRange) {
										dfill[smallerMomIds[j()]].emplace_back(
												momIds[k()]);
									}
								}
							}
						}

						smallerSectionWidth =
								previousSmallerReg.getUpperLimitPoint()[0]
										- previousSmallerReg.getOrigin()[0] - 1;
						// Finish computing the factors
						coefs[2] = coefs[2] / ((double) sectionWidth);
						coefs[3] = coefs[3] / ((double) sectionWidth);
						coefs[4] = coefs[4] / ((double) sectionWidth);
						coefs[5] = coefs[5] / ((double) sectionWidth);
						coefs[6] = coefs[6] / ((double) smallerSectionWidth);
						coefs[7] = coefs[7] / ((double) smallerSectionWidth);
						coefs[8] = coefs[8] / ((double) smallerSectionWidth);
						coefs[9] = coefs[9] / ((double) smallerSectionWidth);
						// Compute the fraction rate
						auto radius = cluster.getReactionRadius();
						double fractionRate = (a1 * exp(-b1 * radius)
								+ (y0 - a1) / (1.0 + c * pow(radius, 2.0))
										* exp(-b2 * pow(radius, 2.0))) * 1.0e-4;
						// Add the size to the vector
						sizeVec.emplace_back(id, previousSmallerId,
								fractionRate, coefs);

						// Reinitialize
						coefs.Init(0.0);
					}

					// Update the previous smaller cluster and the coefs
					previousSmallerId = smallerCluster.getId();
					double distance = (double) i - numXe;
					coefs[0] += 1.0;
					coefs[1] += distance;
					coefs[2] += 1.0;
					coefs[3] += factor;
					coefs[4] += distance;
					coefs[5] += distance * factor;
					coefs[6] += smallerFactor;
					coefs[7] += distance * smallerFactor;
					coefs[8] += 1.0;
					coefs[9] += distance;
				}

				if (previousSmallerId != plsm::invalid<std::size_t>) {
					// Set the connectivities
					dfill[id].emplace_back(id);
					if (clReg.volume() > 1) {
						const auto &momIds = cluster.getMomentIds();
						for (auto j : speciesRange) {
							dfill[id].emplace_back(momIds[j()]);
							dfill[momIds[j()]].emplace_back(id);
							for (auto k : speciesRange) {
								dfill[momIds[j()]].emplace_back(momIds[k()]);
							}
						}
					}
					dfill[xenonId].emplace_back(id);
					if (clReg.volume() > 1) {
						const auto &momIds = cluster.getMomentIds();
						for (auto j : speciesRange) {
							dfill[xenonId].emplace_back(momIds[j()]);
						}
					}
					const auto &previousSmallerReg = neNetwork->getCluster(
							previousSmallerId, plsm::onHost).getRegion();
					dfill[previousSmallerId].emplace_back(id);
					if (clReg.volume() > 1) {
						const auto &momIds = cluster.getMomentIds();
						for (auto j : speciesRange) {
							dfill[previousSmallerId].emplace_back(momIds[j()]);
							if (previousSmallerReg.volume() > 1) {
								const auto &smallerMomIds =
										neNetwork->getCluster(previousSmallerId,
												plsm::onHost).getMomentIds();
								dfill[smallerMomIds[j()]].emplace_back(id);
								for (auto k : speciesRange) {
									dfill[smallerMomIds[j()]].emplace_back(
											momIds[k()]);
								}
							}
						}
					}

					smallerSectionWidth =
							previousSmallerReg.getUpperLimitPoint()[0]
									- previousSmallerReg.getOrigin()[0] - 1;
					// Finish computing the factors
					coefs[2] = coefs[2] / ((double) sectionWidth);
					coefs[3] = coefs[3] / ((double) sectionWidth);
					coefs[4] = coefs[4] / ((double) sectionWidth);
					coefs[5] = coefs[5] / ((double) sectionWidth);
					coefs[6] = coefs[6] / ((double) smallerSectionWidth);
					coefs[7] = coefs[7] / ((double) smallerSectionWidth);
					coefs[8] = coefs[8] / ((double) smallerSectionWidth);
					coefs[9] = coefs[9] / ((double) smallerSectionWidth);
					// Compute the fraction rate
					auto radius = cluster.getReactionRadius();
					double fractionRate = (a1 * exp(-b1 * radius)
							+ (y0 - a1) / (1.0 + c * pow(radius, 2.0))
									* exp(-b2 * pow(radius, 2.0))) * 1.0e-4;
					// Add the size to the vector
					sizeVec.emplace_back(id, previousSmallerId, fractionRate,
							coefs);
				}
			} else {
				// Get the smaller cluster
				comp[NetworkType::Species::Xe] = size - resolutedSize;
				auto smallerCluster = neNetwork->findCluster(comp,
						plsm::onHost);

				// Set the connectivities
				dfill[id].emplace_back(id);
				dfill[xenonId].emplace_back(id);
				dfill[smallerCluster.getId()].emplace_back(id);

				// Add the size to the vector
				coefs[0] = 1.0, coefs[2] = 1.0;
				// Compute the fraction rate
				auto radius = cluster.getReactionRadius();
				double fractionRate = (a1 * exp(-b1 * radius)
						+ (y0 - a1) / (1.0 + c * pow(radius, 2.0))
								* exp(-b2 * pow(radius, 2.0))) * 1.0e-4;
				// Add the size to the vector
				sizeVec.emplace_back(id, smallerCluster.getId(), fractionRate,
						coefs);
			}
		}

	}

	// Print sizeVec
	std::ofstream outputFile;
	outputFile.open("resolutionRateVSradius.txt");
	// Loop on the re-soluting clusters
	for (const auto &currPair : sizeVec) {
		// Get the larger cluster
		auto cluster = neNetwork->getCluster(currPair.larger, plsm::onHost);
		int size = cluster.getRegion().getOrigin()[0];
		double radius = cluster.getReactionRadius();

		outputFile << size << " " << radius << " " << currPair.fractionRate
				<< std::endl;
	}
	outputFile.close();

	return;
}

void ReSolutionHandler::updateReSolutionRate(double rate) {
	// Scale the fission rate with the standard one
	resolutionRate = 1.0e8 * rate / fissionYield;

	return;
}

void ReSolutionHandler::setFissionYield(double yield) {
	fissionYield = yield;

	return;
}

void ReSolutionHandler::computeReSolution(
		experimental::IReactionNetwork &network, double *concOffset,
		double *updatedConcOffset, int xi, int xs, int yj, int zk) {

	using NetworkType =
	experimental::NEReactionNetwork;
	auto neNetwork = dynamic_cast<NetworkType*>(&network);
	// Get the single xenon which connects to every cluster
	NetworkType::Composition comp;
	for (auto i : neNetwork->getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[NetworkType::Species::Xe] = 1;
	auto singleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto xenonId = singleXenon.getId();

	// Loop on the re-soluting clusters
	for (const auto &currPair : sizeVec) {
		// Get the larger cluster
		auto cluster = neNetwork->getCluster(currPair.larger, plsm::onHost);
		auto momId = cluster.getMomentIds()[0];
		double rate = currPair.fractionRate * resolutionRate;
		// Get the re-solution cluster
		auto resoCluster = neNetwork->getCluster(currPair.smaller,
				plsm::onHost);
		auto resoMomId = resoCluster.getMomentIds()[0];

		// Get the initial concentration of the larger xenon cluster
		double l0 = concOffset[currPair.larger], l1 = concOffset[momId];
		// Update the concentrations
		updatedConcOffset[currPair.larger] -= rate
				* (currPair.coefs[2] * l0 + currPair.coefs[4] * l1);
		updatedConcOffset[momId] -= rate
				* (currPair.coefs[3] * l0 + currPair.coefs[5] * l1);
		if (resoCluster.getRegion().volume() > 1) {
			updatedConcOffset[currPair.smaller] += rate
					* (currPair.coefs[8] * l0 + currPair.coefs[9] * l1);
			updatedConcOffset[resoMomId] += rate
					* (currPair.coefs[6] * l0 + currPair.coefs[7] * l1);
		} else {
			updatedConcOffset[currPair.smaller] += rate
					* (currPair.coefs[0] * l0 + currPair.coefs[1] * l1);
		}
		updatedConcOffset[xenonId] += rate
				* (currPair.coefs[0] * l0 + currPair.coefs[1] * l1);
	}

	return;
}

int ReSolutionHandler::computePartialsForReSolution(
		experimental::IReactionNetwork &network, double *val, int *indices,
		int xi, int xs, int yj, int zk) {

	using NetworkType =
	experimental::NEReactionNetwork;
	auto neNetwork = dynamic_cast<NetworkType*>(&network);
	// Get the single xenon which connects to every cluster
	NetworkType::Composition comp;
	for (auto i : neNetwork->getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[NetworkType::Species::Xe] = 1;
	auto singleXenon = neNetwork->findCluster(comp, plsm::onHost);
	auto xenonId = singleXenon.getId();

	// Loop on the re-soluting clusters
	int i = 0;
	for (const auto &currPair : sizeVec) {
		// Get the larger cluster
		auto cluster = neNetwork->getCluster(currPair.larger, plsm::onHost);
		auto momId = cluster.getMomentIds()[0];
		double rate = currPair.fractionRate * resolutionRate;
		// Get the re-solution cluster
		auto resoCluster = neNetwork->getCluster(currPair.smaller,
				plsm::onHost);
		auto resoMomId = resoCluster.getMomentIds()[0];

		// Set the partial derivatives
		auto baseIndex = i * 10;
		indices[5 * i] = currPair.larger;
		val[baseIndex] = -rate * currPair.coefs[2];
		indices[(5 * i) + 1] = momId;
		val[(baseIndex) + 1] = -rate * currPair.coefs[4];
		val[(baseIndex) + 2] = -rate * currPair.coefs[3];
		val[(baseIndex) + 3] = -rate * currPair.coefs[5]; // Large cluster
		if (resoCluster.getRegion().volume() > 1) {
			indices[(5 * i) + 2] = currPair.smaller;
			val[(baseIndex) + 4] = rate * currPair.coefs[8];
			indices[(5 * i) + 3] = resoMomId;
			val[(baseIndex) + 5] = rate * currPair.coefs[9];
			val[(baseIndex) + 6] = rate * currPair.coefs[6];
			val[(baseIndex) + 7] = rate * currPair.coefs[7]; // Smaller cluster
		} else {
			indices[(5 * i) + 2] = currPair.smaller;
			val[(baseIndex) + 4] = rate * currPair.coefs[0];
			indices[(5 * i) + 3] = resoMomId;
			val[(baseIndex) + 5] = rate * currPair.coefs[1];
			val[(baseIndex) + 6] = 0.0;
			val[(baseIndex) + 7] = 0.0; // Smaller cluster
		}
		indices[(5 * i) + 4] = xenonId;
		val[(baseIndex) + 8] = rate * currPair.coefs[0];
		val[(baseIndex) + 9] = rate * currPair.coefs[1]; // Xe_1

		// Increment i
		i++;
	}

	return sizeVec.size();
}

}/* end namespace xolotlCore */

