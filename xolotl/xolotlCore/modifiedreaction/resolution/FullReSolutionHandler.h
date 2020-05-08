#ifndef FULLRESOLUTIONHANDLER_H
#define FULLRESOLUTIONHANDLER_H

// Includes
#include <ReSolutionHandler.h>
#include <Constants.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <NESuperCluster.h>

namespace xolotlCore {

/**
 * This class realizes IReSolutionHandler interface, responsible for the
 * re-solution of xenon clusters with full destruction.
 */
class FullReSolutionHandler: public ReSolutionHandler {

public:

	/**
	 * The constructor
	 */
	FullReSolutionHandler() :
			ReSolutionHandler() {
	}

	/**
	 * The destructor
	 */
	~FullReSolutionHandler() {
	}

	/**
	 * The initialize method has to add connectivity between the Xe clusters,
	 * It must also initialize the rates of the reactions.
	 *
	 * \see IReSolutionHandler.h
	 */
	void initialize(const IReactionNetwork &network,
			double electronicStoppingPower) {
		// Get the single xenon which connects to every cluster
		auto singleXenon = network.get(Species::Xe, 1);

		if (!singleXenon) {
			// Inform the user
			std::cout << "The re-solution won't happen because "
					"the single xenon cluster is missing." << std::endl;

			return;
		}

		const double lambda = 1.0;
		const double mu = 6000.0;

		// Loop on the clusters
		auto allClusters = network.getAll();
		std::for_each(allClusters.begin(), allClusters.end(),
				[&network, this, &singleXenon, &lambda, &mu](
						IReactant &cluster) {
					// Get its size
					auto size = cluster.getSize();
					if (size > 1) {

						// Initialize the coefficients
						Array<double, 10> coefs;
						coefs.Init(0.0);
						// Different if it is a super cluster
						if (cluster.getType() == ReactantType::NESuper) {
							// Loop the bounds of the super cluster
							auto &superCluster =
									static_cast<NESuperCluster&>(cluster);
							int sectionWidth = superCluster.getSectionWidth();
							double dispersion = superCluster.getDispersion();
							double numXe = superCluster.getAverage();
							for (int i = numXe + 1.0
									- (double) sectionWidth / 2.0;
									i <= numXe + (double) sectionWidth / 2.0;
									i++) {
								// Update the quantities needed for the coefs
								double factor = ((double) i - numXe)
										/ dispersion, distance =
										superCluster.getDistance(i);

								// Update the coefs
								coefs[0] += 1.0;
								coefs[1] += distance;
								coefs[2] += 1.0;
								coefs[3] += factor;
								coefs[4] += distance;
								coefs[5] += distance * factor;
								coefs[6] += (double) i;
							}

							// Set the connectivities
							cluster.setDissociationConnectivity(
									cluster.getId());
							cluster.setDissociationConnectivity(
									cluster.getMomentId());
							singleXenon->setDissociationConnectivity(
									cluster.getId());
							singleXenon->setDissociationConnectivity(
									cluster.getMomentId());
							// Finish computing the factors
							coefs[2] = coefs[2] / ((double) sectionWidth);
							coefs[3] = coefs[3] / ((double) sectionWidth);
							coefs[4] = coefs[4] / ((double) sectionWidth);
							coefs[5] = coefs[5] / ((double) sectionWidth);
							coefs[6] = coefs[6] / ((double) sectionWidth);
							// Compute the fraction rate
							auto radius = cluster.getReactionRadius();
							double fractionRate = 4.0 * xolotlCore::pi * mu
									* pow(radius + lambda, 2.0);
							// Add the size to the vector
							sizeVec.emplace_back(&cluster, nullptr,
									fractionRate, coefs);
						} else {
							// Set the connectivities
							cluster.setDissociationConnectivity(
									cluster.getId());
							singleXenon->setDissociationConnectivity(
									cluster.getId());
							// Add the size to the vector
							coefs[0] = 1.0, coefs[2] = 1.0, coefs[6] = size;
							// Compute the fraction rate
							auto radius = cluster.getReactionRadius();
							double fractionRate = 4.0 * xolotlCore::pi * mu
									* pow(radius + lambda, 2.0);
							// Add the size to the vector
							sizeVec.emplace_back(&cluster, nullptr,
									fractionRate, coefs);
						}
					}
				});

		// Print sizeVec
		std::ofstream outputFile;
		outputFile.open("resolutionRateVSradius.txt");
		// Loop on the re-soluting clusters
		for (const auto &currPair : sizeVec) {
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

	void updateReSolutionRate(double rate) {
		// Simply save the fission rate
		resolutionRate = rate / fissionYield;

		return;
	}

	/**
	 * Compute the flux due to the re-solution for all the cluster,
	 * given the position index xi.
	 * This method is called by the RHSFunction from the PetscSolver.
	 *
	 * Xe_i --> i Xe_1
	 *
	 * F(Xe_i) = -1/i * F(Xe_1) = -rate * C_(Xe_i)
	 *
	 * \see IReSolutionHandler.h
	 */
	void computeReSolution(const IReactionNetwork &network, double *concOffset,
			double *updatedConcOffset, int xi, int xs, int yj = 0, int zk = 0) {
		// Get the single xenon which connects to every cluster
		auto singleXenon = network.get(Species::Xe, 1);
		int xenonId = singleXenon->getId() - 1;

		// Loop on the re-soluting clusters
		for (const auto &currPair : sizeVec) {
			// Get the larger cluster
			auto cluster = currPair.larger;
			int id = cluster->getId() - 1;
			int momId = cluster->getMomentId() - 1;
			double size = currPair.coefs[6];
			double rate = currPair.fractionRate * resolutionRate;

			// Get the initial concentration of the larger xenon cluster
			double l0 = concOffset[id], l1 = concOffset[momId];
			// Update the concentrations
			updatedConcOffset[id] -= rate
					* (currPair.coefs[2] * l0 + currPair.coefs[4] * l1);
			updatedConcOffset[momId] -= rate
					* (currPair.coefs[3] * l0 + currPair.coefs[5] * l1);
			updatedConcOffset[xenonId] += size * rate
					* (currPair.coefs[0] * l0 + currPair.coefs[1] * l1);
		}

		return;
	}

	/**
	 * Compute the partials due to the re-solution for all the
	 * clusters given the position index xi. Returns the number of xenon
	 * clusters that are re-soluted at this grid point.
	 * This method is called by the RHSJacobian from the PetscSolver.
	 *
	 * Xe_i --> i Xe_1
	 *
	 * dF(Xe_i)/dC_(Xe_i) = -1/n * dF(Xe_1)/dC_(Xe_i)
	 * 		= -rate
	 *
	 * \see IReSolutionHandler.h
	 */
	int computePartialsForReSolution(const IReactionNetwork &network,
			double *val, int *indices, int xi, int xs, int yj = 0, int zk = 0) {
		// Get the single xenon which connects to every cluster
		auto singleXenon = network.get(Species::Xe, 1);
		int xenonId = singleXenon->getId() - 1;

		// Loop on the re-soluting clusters
		int i = 0;
		for (const auto &currPair : sizeVec) {
			// Get the larger cluster
			auto cluster = currPair.larger;
			int id = cluster->getId() - 1;
			int momId = cluster->getMomentId() - 1;
			double size = currPair.coefs[6];
			double rate = currPair.fractionRate * resolutionRate;

			// Set the partial derivatives
			auto baseIndex = i * 10;
			indices[5 * i] = id;
			val[baseIndex] = -rate * currPair.coefs[2];
			indices[(5 * i) + 1] = momId;
			val[(baseIndex) + 1] = -rate * currPair.coefs[4];
			val[(baseIndex) + 2] = -rate * currPair.coefs[3];
			val[(baseIndex) + 3] = -rate * currPair.coefs[5]; // Large cluster

			indices[(5 * i) + 2] = 0;
			val[(baseIndex) + 4] = 0.0;
			indices[(5 * i) + 3] = 0;
			val[(baseIndex) + 5] = 0.0;
			val[(baseIndex) + 6] = 0.0;
			val[(baseIndex) + 7] = 0.0; // Smaller cluster

			indices[(5 * i) + 4] = xenonId;
			val[(baseIndex) + 8] = size * rate * currPair.coefs[0];
			val[(baseIndex) + 9] = size * rate * currPair.coefs[1]; // Xe_1

			// Increment i
			i++;
		}

		return sizeVec.size();
	}
};
//end class FullReSolutionHandler

} /* namespace xolotlCore */
#endif
