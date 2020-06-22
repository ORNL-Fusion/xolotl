#ifndef PARTIALRESOLUTIONHANDLER_H
#define PARTIALRESOLUTIONHANDLER_H

// Includes
#include <Constants.h>
#include <NESuperCluster.h>
#include <ReSolutionHandler.h>

#include <algorithm>
#include <fstream>
#include <iostream>

namespace xolotlCore
{
/**
 * This class realizes IReSolutionHandler interface, responsible for the
 * re-solution of xenon clusters on a partial basis.
 */
class PartialReSolutionHandler : public ReSolutionHandler
{
public:
	/**
	 * The constructor
	 */
	PartialReSolutionHandler() : ReSolutionHandler()
	{
	}

	/**
	 * The destructor
	 */
	~PartialReSolutionHandler()
	{
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reactions.
	 *
	 * \see IReSolutionHandler.h
	 */
	void
	initialize(const IReactionNetwork& network, double electronicStoppingPower)
	{
		// Get the single xenon which connects to every cluster
		auto singleXenon = network.get(Species::Xe, 1);

		if (!singleXenon) {
			// Inform the user
			std::cout << "The re-solution won't happen because "
						 "the single xenon cluster is missing."
					  << std::endl;

			return;
		}

		// Set the fit variables depending on the electronic stopping power
		double y0 = 0.0, a1 = 0.0, b1 = 0.0, b2 = 0.0, c = 0.0;
		if (electronicStoppingPower > 0.87) {
			y0 = 11.0851, a1 = 1.5052, b1 = 0.0362, b2 = 0.0203, c = 3.4123;
		}
		else if (electronicStoppingPower > 0.82) {
			y0 = 10.6297, a1 = 1.3479, b1 = 0.0438, b2 = 0.0241, c = 4.2214;
		}
		else if (electronicStoppingPower > 0.77) {
			y0 = 10.1521, a1 = 1.1986, b1 = 0.0546, b2 = 0.0299, c = 5.4612;
		}
		else if (electronicStoppingPower > 0.71) {
			y0 = 9.1816, a1 = 0.949, b1 = 0.0703, b2 = 0.0371, c = 7.982;
		}
		else if (electronicStoppingPower > 0.67) {
			y0 = 8.6745, a1 = 0.8401, b1 = 0.0792, b2 = 0.0407, c = 9.6585;
		}
		else if (electronicStoppingPower > 0.62) {
			y0 = 7.6984, a1 = 0.6721, b1 = 0.1028, b2 = 0.0526, c = 14.272;
		}
		else if (electronicStoppingPower > 0.57) {
			y0 = 6.3925, a1 = 0.5025, b1 = 0.1411, b2 = 0.0727, c = 23.1967;
		}
		else if (electronicStoppingPower > 0.52) {
			y0 = 4.6175, a1 = 0.3433, b1 = 0.2284, b2 = 0.1276, c = 45.6624;
		}
		else {
			y0 = 2.3061, a1 = 0.2786, b1 = 1.1008, b2 = 1.605, c = 150.6689;
		}

		// Compute the thermal spike energy
		const double seeff = 16.0;

		// Loop on the clusters
		auto allClusters = network.getAll();
		std::for_each(allClusters.begin(), allClusters.end(),
			[&network, this, &singleXenon, &y0, &a1, &b1, &b2, &c, &seeff](
				IReactant& cluster) {
				// Get its size
				auto size = cluster.getSize();
				auto radius = cluster.getReactionRadius();

				// The re-soluted size is always 1
				int resolutedSize = (double)size * 0.25 *
						(1.0 -
							exp(-0.05 * (seeff - 9.04) / pow(radius, 1.47))) +
					0.5;
				resolutedSize = std::max(1, resolutedSize);

				// Partial case
				if (resolutedSize > 0 && resolutedSize < size - 1 && size > 1) {
					// Initialize the coefficients
					Array<double, 11> coefs;
					coefs.Init(0.0);
					// Different if it is a super cluster
					if (cluster.getType() == ReactantType::NESuper) {
						// Loop the bounds of the super cluster
						auto& superCluster =
							static_cast<NESuperCluster&>(cluster);
						int sectionWidth = superCluster.getSectionWidth(),
							smallerSectionWidth = 0.0;
						double dispersion = superCluster.getDispersion();
						double numXe = superCluster.getAverage();
						IReactant* previousSmaller = nullptr;
						for (int i = numXe + 1.0 - (double)sectionWidth / 2.0;
							 i <= numXe + (double)sectionWidth / 2.0; i++) {
							// Update the quantities needed for the coefs
							double factor = ((double)i - numXe) / dispersion,
								   distance = superCluster.getDistance(i),
								   smallerFactor = 0.0;

							// Get the smaller cluster
							auto smallerCluster =
								network.get(Species::Xe, i - resolutedSize);
							// If it exists
							if (!smallerCluster) {
								// Get the super cluster
								smallerCluster = network.getSuperFromComp(
									i - resolutedSize, 0, 0, 0);
								auto smallerSuper =
									static_cast<NESuperCluster*>(
										smallerCluster);
								// Update the corresponding coefs
								smallerFactor =
									((double)(i - resolutedSize) -
										smallerSuper->getAverage()) /
									smallerSuper->getDispersion();
							}

							// Check if it was the same as before
							if (smallerCluster != previousSmaller &&
								previousSmaller) {
								// Set the connectivities
								cluster.setDissociationConnectivity(
									cluster.getId());
								cluster.setDissociationConnectivity(
									cluster.getMomentId());
								singleXenon->setDissociationConnectivity(
									cluster.getId());
								singleXenon->setDissociationConnectivity(
									cluster.getMomentId());
								previousSmaller->setDissociationConnectivity(
									cluster.getId());
								previousSmaller->setDissociationConnectivity(
									cluster.getMomentId());
								smallerSectionWidth =
									((NECluster*)previousSmaller)
										->getSectionWidth();
								// Finish computing the factors
								coefs[2] = coefs[2] / ((double)sectionWidth);
								coefs[3] = coefs[3] / ((double)sectionWidth);
								coefs[4] = coefs[4] / ((double)sectionWidth);
								coefs[5] = coefs[5] / ((double)sectionWidth);
								coefs[6] =
									coefs[6] / ((double)smallerSectionWidth);
								coefs[7] =
									coefs[7] / ((double)smallerSectionWidth);
								coefs[8] =
									coefs[8] / ((double)smallerSectionWidth);
								coefs[9] =
									coefs[9] / ((double)smallerSectionWidth);
								coefs[10] = coefs[10] / ((double)coefs[0]);
								// Compute the fraction rate
								double fractionRate =
									(a1 * exp(-b1 * radius) +
										(y0 - a1) /
											(1.0 + c * pow(radius, 2.0)) *
											exp(-b2 * pow(radius, 2.0))) *
									1.0e-4 / (double)resolutedSize;
								// Add the size to the vector
								sizeVec.emplace_back(&cluster, previousSmaller,
									fractionRate, coefs);

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
							coefs[8] += 1.0;
							coefs[9] += distance;
							coefs[10] += (double)resolutedSize;
						}
						if (previousSmaller) {
							// Set the connectivities
							cluster.setDissociationConnectivity(
								cluster.getId());
							cluster.setDissociationConnectivity(
								cluster.getMomentId());
							singleXenon->setDissociationConnectivity(
								cluster.getId());
							singleXenon->setDissociationConnectivity(
								cluster.getMomentId());
							previousSmaller->setDissociationConnectivity(
								cluster.getId());
							previousSmaller->setDissociationConnectivity(
								cluster.getMomentId());
							smallerSectionWidth = ((NECluster*)previousSmaller)
													  ->getSectionWidth();
							// Finish computing the factors
							coefs[2] = coefs[2] / ((double)sectionWidth);
							coefs[3] = coefs[3] / ((double)sectionWidth);
							coefs[4] = coefs[4] / ((double)sectionWidth);
							coefs[5] = coefs[5] / ((double)sectionWidth);
							coefs[6] = coefs[6] / ((double)smallerSectionWidth);
							coefs[7] = coefs[7] / ((double)smallerSectionWidth);
							coefs[8] = coefs[8] / ((double)smallerSectionWidth);
							coefs[9] = coefs[9] / ((double)smallerSectionWidth);
							coefs[10] = coefs[10] / ((double)coefs[0]);
							// Compute the fraction rate
							double fractionRate =
								(a1 * exp(-b1 * radius) +
									(y0 - a1) / (1.0 + c * pow(radius, 2.0)) *
										exp(-b2 * pow(radius, 2.0))) *
								1.0e-4 / (double)resolutedSize;
							// Add the size to the vector
							sizeVec.emplace_back(
								&cluster, previousSmaller, fractionRate, coefs);
						}
					}
					else {
						// Get the smaller cluster
						auto smallerCluster =
							network.get(Species::Xe, size - resolutedSize);

						// Set the connectivities
						cluster.setDissociationConnectivity(cluster.getId());
						singleXenon->setDissociationConnectivity(
							cluster.getId());
						smallerCluster->setDissociationConnectivity(
							cluster.getId());
						// Add the size to the vector
						coefs[0] = 1.0, coefs[2] = 1.0,
						coefs[10] = resolutedSize;
						// Compute the fraction rate
						double fractionRate =
							(a1 * exp(-b1 * radius) +
								(y0 - a1) / (1.0 + c * pow(radius, 2.0)) *
									exp(-b2 * pow(radius, 2.0))) *
							1.0e-4 / (double)resolutedSize;
						// Add the size to the vector
						sizeVec.emplace_back(
							&cluster, smallerCluster, fractionRate, coefs);
					}
				}

				// Full case
				else if (resolutedSize > 0 && size > 1) {
					// Initialize the coefficients
					Array<double, 11> coefs;
					coefs.Init(0.0);
					// Different if it is a super cluster
					if (cluster.getType() == ReactantType::NESuper) {
						// Loop the bounds of the super cluster
						auto& superCluster =
							static_cast<NESuperCluster&>(cluster);
						int sectionWidth = superCluster.getSectionWidth();
						double dispersion = superCluster.getDispersion();
						double numXe = superCluster.getAverage();
						for (int i = numXe + 1.0 - (double)sectionWidth / 2.0;
							 i <= numXe + (double)sectionWidth / 2.0; i++) {
							// Update the quantities needed for the coefs
							double factor = ((double)i - numXe) / dispersion,
								   distance = superCluster.getDistance(i);

							// Update the coefs
							coefs[0] += 1.0;
							coefs[1] += distance;
							coefs[2] += 1.0;
							coefs[3] += factor;
							coefs[4] += distance;
							coefs[5] += distance * factor;
							coefs[10] += (double)i;
						}

						// Set the connectivities
						cluster.setDissociationConnectivity(cluster.getId());
						cluster.setDissociationConnectivity(
							cluster.getMomentId());
						singleXenon->setDissociationConnectivity(
							cluster.getId());
						singleXenon->setDissociationConnectivity(
							cluster.getMomentId());
						// Finish computing the factors
						coefs[2] = coefs[2] / ((double)sectionWidth);
						coefs[3] = coefs[3] / ((double)sectionWidth);
						coefs[4] = coefs[4] / ((double)sectionWidth);
						coefs[5] = coefs[5] / ((double)sectionWidth);
						coefs[10] = coefs[10] / ((double)sectionWidth);
						// Compute the fraction rate
						double fractionRate =
							(a1 * exp(-b1 * radius) +
								(y0 - a1) / (1.0 + c * pow(radius, 2.0)) *
									exp(-b2 * pow(radius, 2.0))) *
							1.0e-4 / (double)size;
						// Add the size to the vector
						sizeVec.emplace_back(
							&cluster, nullptr, fractionRate, coefs);
					}
					else {
						// Set the connectivities
						cluster.setDissociationConnectivity(cluster.getId());
						singleXenon->setDissociationConnectivity(
							cluster.getId());
						// Add the size to the vector
						coefs[0] = 1.0, coefs[2] = 1.0, coefs[10] = size;
						// Compute the fraction rate
						double fractionRate =
							(a1 * exp(-b1 * radius) +
								(y0 - a1) / (1.0 + c * pow(radius, 2.0)) *
									exp(-b2 * pow(radius, 2.0))) *
							1.0e-4 / (double)size;
						// Add the size to the vector
						sizeVec.emplace_back(
							&cluster, nullptr, fractionRate, coefs);
					}
				}
			});

		// Print sizeVec
		std::ofstream outputFile;
		outputFile.open("resolutionRateVSradius.txt");
		// Loop on the re-soluting clusters
		for (const auto& currPair : sizeVec) {
			// Get the larger cluster
			auto cluster = currPair.larger;
			int size = cluster->getSize();
			double radius = cluster->getReactionRadius();

			outputFile << size << " " << radius << " " << currPair.fractionRate
					   << std::endl;
		}
		outputFile.close();

		return;
	}

	void
	updateReSolutionRate(double rate)
	{
		// Scale the fission rate with the standard one
		resolutionRate = 1.0e8 * rate / fissionYield;

		return;
	}

	/**
	 * Compute the flux due to the re-solution for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * Xe_i --> Xe_(i-n) + n Xe_1 (n=1)
	 *
	 * F(Xe_i) = -F[Xe_(i-n)] = -1/n * F(Xe_1) = -rate * C_(Xe_i)
	 *
	 * \see IReSolutionHandler.h
	 */
	void
	computeReSolution(const IReactionNetwork& network, double* concOffset,
		double* updatedConcOffset, int xi, int xs, int yj = 0, int zk = 0)
	{
		// Get the single xenon which connects to every cluster
		auto singleXenon = network.get(Species::Xe, 1);
		int xenonId = singleXenon->getId() - 1;

		// Loop on the re-soluting clusters
		for (const auto& currPair : sizeVec) {
			// Get the larger cluster
			auto cluster = currPair.larger;
			int id = cluster->getId() - 1;
			int momId = cluster->getMomentId() - 1;
			double size = currPair.coefs[10];
			double rate = currPair.fractionRate * resolutionRate;

			// Get the initial concentration of the larger xenon cluster
			double l0 = concOffset[id], l1 = concOffset[momId];
			// Update the concentrations
			updatedConcOffset[id] -=
				rate * (currPair.coefs[2] * l0 + currPair.coefs[4] * l1);
			updatedConcOffset[momId] -=
				rate * (currPair.coefs[3] * l0 + currPair.coefs[5] * l1);
			// Get the re-solution cluster
			if (currPair.smaller) {
				auto resoCluster = currPair.smaller;
				int resoId = resoCluster->getId() - 1;
				int resoMomId = resoCluster->getMomentId() - 1;
				if (resoCluster->getType() == ReactantType::NESuper) {
					updatedConcOffset[resoId] += rate *
						(currPair.coefs[8] * l0 + currPair.coefs[9] * l1);
					updatedConcOffset[resoMomId] += rate *
						(currPair.coefs[6] * l0 + currPair.coefs[7] * l1);
				}
				else {
					updatedConcOffset[resoId] += rate *
						(currPair.coefs[0] * l0 + currPair.coefs[1] * l1);
				}
			}
			updatedConcOffset[xenonId] +=
				size * rate * (currPair.coefs[0] * l0 + currPair.coefs[1] * l1);
		}

		return;
	}

	/**
	 * Compute the partials due to the re-solution for all the
	 * clusters given the position index xi. Returns the number of xenon
	 * clusters that are re-soluted at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Xe_i --> Xe_(i-n) + n Xe_1 (n=1)
	 *
	 * dF(Xe_i)/dC_(Xe_i) = -dF[Xe_(i-n)]/dC_(Xe_i) = -1/n * dF(Xe_1)/dC_(Xe_i)
	 * 		= -rate
	 *
	 * \see IReSolutionHandler.h
	 */
	int
	computePartialsForReSolution(const IReactionNetwork& network, double* val,
		int* indices, int xi, int xs, int yj = 0, int zk = 0)
	{
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
			double size = currPair.coefs[10];
			double rate = currPair.fractionRate * resolutionRate;

			// Set the partial derivatives
			auto baseIndex = i * 10;
			indices[5 * i] = id;
			val[baseIndex] = -rate * currPair.coefs[2];
			indices[(5 * i) + 1] = momId;
			val[(baseIndex) + 1] = -rate * currPair.coefs[4];
			val[(baseIndex) + 2] = -rate * currPair.coefs[3];
			val[(baseIndex) + 3] = -rate * currPair.coefs[5]; // Large cluster
			// Get the re-solution cluster
			if (currPair.smaller) {
				auto resoCluster = currPair.smaller;
				int resoId = resoCluster->getId() - 1;
				int resoMomId = resoCluster->getMomentId() - 1;

				if (resoCluster->getType() == ReactantType::NESuper) {
					indices[(5 * i) + 2] = resoId;
					val[(baseIndex) + 4] = rate * currPair.coefs[8];
					indices[(5 * i) + 3] = resoMomId;
					val[(baseIndex) + 5] = rate * currPair.coefs[9];
					val[(baseIndex) + 6] = rate * currPair.coefs[6];
					val[(baseIndex) + 7] =
						rate * currPair.coefs[7]; // Smaller cluster
				}
				else {
					indices[(5 * i) + 2] = resoId;
					val[(baseIndex) + 4] = rate * currPair.coefs[0];
					indices[(5 * i) + 3] = resoMomId;
					val[(baseIndex) + 5] = rate * currPair.coefs[1];
					val[(baseIndex) + 6] = 0.0;
					val[(baseIndex) + 7] = 0.0; // Smaller cluster
				}
			}
			else {
				indices[(5 * i) + 2] = 0;
				val[(baseIndex) + 4] = 0.0;
				indices[(5 * i) + 3] = 0;
				val[(baseIndex) + 5] = 0.0;
				val[(baseIndex) + 6] = 0.0;
				val[(baseIndex) + 7] = 0.0; // Smaller cluster
			}
			indices[(5 * i) + 4] = xenonId;
			val[(baseIndex) + 8] = size * rate * currPair.coefs[0];
			val[(baseIndex) + 9] = size * rate * currPair.coefs[1]; // Xe_1

			// Increment i
			i++;
		}

		return sizeVec.size();
	}
};
// end class PartialReSolutionHandler

} /* namespace xolotlCore */
#endif
