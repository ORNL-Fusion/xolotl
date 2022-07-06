// Includes
#include <array>
#include <iostream>

#include <xolotl/core/Constants.h>
#include <xolotl/core/modified/SoretDiffusionHandler.h>

namespace xolotl
{
namespace core
{
namespace modified
{
void
SoretDiffusionHandler::computeDiffusion(network::IReactionNetwork& network,
	double** concVector, double* updatedConcOffset, double hxLeft,
	double hxRight, int ix, double, int, double, int) const
{
	int diffClusterIdx = 0;

	for (auto const& currId : diffusingClusters) {
		auto cluster = network.getClusterCommon(currId);
		// Get the initial concentrations
		double oldConc = concVector[0][currId];
		double oldLeftConc = concVector[1][currId];
		double oldRightConc = concVector[2][currId];

		double leftTemp = cluster.getTemperature(ix),
			   midTemp = cluster.getTemperature(ix + 1),
			   rightTemp = cluster.getTemperature(ix + 2);
		double leftDiff = cluster.getDiffusionCoefficient(ix),
			   midDiff = cluster.getDiffusionCoefficient(ix + 1),
			   rightDiff = cluster.getDiffusionCoefficient(ix + 2);

		updatedConcOffset[currId] -= 2.0 * beta[diffClusterIdx] * midDiff *
				oldConc *
				(leftTemp + (hxLeft / hxRight) * rightTemp -
					(1.0 + (hxLeft / hxRight)) * midTemp) /
				(hxLeft * (hxLeft + hxRight)) +
			beta[diffClusterIdx] * midDiff * (oldRightConc - oldLeftConc) *
				(rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)) +
			beta[diffClusterIdx] * oldConc * (rightDiff - leftDiff) *
				(rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight));

		// Increase the index
		diffClusterIdx++;
	}

	return;
}

bool
SoretDiffusionHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, double** concVector, double* val,
	IdType* indices, double hxLeft, double hxRight, int ix, double, int, double,
	int) const
{
	int diffClusterIdx = 0;

	for (auto const& currId : diffusingClusters) {
		auto cluster = network.getClusterCommon(currId);
		// Set the cluster index, the PetscSolver will use it to compute
		// the row and column indices for the Jacobian
		indices[diffClusterIdx] = currId;
		// Get the initial concentrations
		double oldConc = concVector[0][currId];
		double oldLeftConc = concVector[1][currId];
		double oldRightConc = concVector[2][currId];

		double leftTemp = cluster.getTemperature(ix),
			   midTemp = cluster.getTemperature(ix + 1),
			   rightTemp = cluster.getTemperature(ix + 2);
		double leftDiff = cluster.getDiffusionCoefficient(ix),
			   midDiff = cluster.getDiffusionCoefficient(ix + 1),
			   rightDiff = cluster.getDiffusionCoefficient(ix + 2);

		// Compute the partial derivatives for diffusion of this cluster
		// for the middle, left, and right grid point
		val[diffClusterIdx * 6] = -2.0 * beta[diffClusterIdx] * midDiff *
				(leftTemp + (hxLeft / hxRight) * rightTemp -
					(1.0 + (hxLeft / hxRight)) * midTemp) /
				(hxLeft * (hxLeft + hxRight)) -
			beta[diffClusterIdx] * (rightDiff - leftDiff) *
				(rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)); // middle conc

		val[(diffClusterIdx * 6) + 1] = beta[diffClusterIdx] * midDiff *
			(rightTemp - leftTemp) /
			((hxLeft + hxRight) * (hxLeft + hxRight)); // left conc

		val[(diffClusterIdx * 6) + 2] = -beta[diffClusterIdx] * midDiff *
			(rightTemp - leftTemp) /
			((hxLeft + hxRight) * (hxLeft + hxRight)); // right conc

		val[(diffClusterIdx * 6) + 3] = 0.0; // middle temp

		val[(diffClusterIdx * 6) + 4] = 0.0; // left temp

		val[(diffClusterIdx * 6) + 5] = 0.0; // right temp

		// Increase the index
		diffClusterIdx++;
	}

	return true;
}

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
