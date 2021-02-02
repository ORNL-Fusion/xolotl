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
double
SoretDiffusionHandler::getLocalHeatFactor(int xi) const
{
	//	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	//	if (x < 100.0)
	//		return 0.2;
	//	//	if (x < 125.0) return 0.016 * x - 1.0;
	return 1.0;
}

void
SoretDiffusionHandler::computeDiffusion(network::IReactionNetwork& network,
	double** concVector, double* updatedConcOffset, double hxLeft,
	double hxRight, int ix, double, int, double, int) const
{
	// Adjust the constants to the material
	double localHeatCond = getLocalHeatFactor(ix) * heatConductivity;

	if (ix == surfacePosition) {
		for (auto const& currId : diffusingClusters) {
			auto cluster = network.getClusterCommon(currId);
			// Get the initial concentrations
			double oldConc = concVector[0][currId];
			double oldLeftConc = concVector[1][currId];
			double oldRightConc = concVector[2][currId];
			double midTemp = concVector[0][dof], rightTemp = concVector[2][dof];
			double midDiff = cluster.getDiffusionCoefficient(ix + 1),
				   rightDiff = cluster.getDiffusionCoefficient(ix + 2),
				   leftDiff = cluster.getDiffusionFactor() *
				exp(-cluster.getMigrationEnergy() /
					(kBoltzmann *
						(rightTemp +
							(hxLeft + hxRight) * heatFlux / localHeatCond)));

			updatedConcOffset[currId] += 2.0 *
					(J +
						(beta * midDiff * oldConc * heatFlux) / localHeatCond) /
					hxLeft +
				2.0 * midDiff * (oldRightConc - oldConc) / (hxLeft * hxRight) -
				(J + (beta * midDiff * oldConc * heatFlux) / localHeatCond) *
					(rightDiff - leftDiff) / (midDiff * (hxLeft + hxRight));

			// Second part
			updatedConcOffset[currId] -= 2.0 * beta * midDiff * oldConc *
					heatFlux / (hxLeft * localHeatCond) +
				2.0 * beta * midDiff * oldConc * (rightTemp - midTemp) /
					(hxLeft * hxRight) +
				beta * heatFlux *
					(J +
						(beta * midDiff * oldConc * heatFlux) / localHeatCond) /
					localHeatCond -
				beta * oldConc * heatFlux * (rightDiff - leftDiff) /
					(localHeatCond * (hxLeft + hxRight));
		}
	}
	else {
		for (auto const& currId : diffusingClusters) {
			auto cluster = network.getClusterCommon(currId);
			// Get the initial concentrations
			double oldConc = concVector[0][currId];
			double oldLeftConc = concVector[1][currId];
			double oldRightConc = concVector[2][currId];

			double J = 1.0e3; // nm-2 s-1
			//			double J = 0.0; // nm-2 s-1
			double leftTemp = concVector[1][dof], midTemp = concVector[0][dof],
				   rightTemp = concVector[2][dof];
			double leftDiff = cluster.getDiffusionCoefficient(ix),
				   midDiff = cluster.getDiffusionCoefficient(ix + 1),
				   rightDiff = cluster.getDiffusionCoefficient(ix + 2);

			updatedConcOffset[currId] -= 2.0 * beta * midDiff * oldConc *
					(leftTemp + (hxLeft / hxRight) * rightTemp -
						(1.0 + (hxLeft / hxRight)) * midTemp) /
					(hxLeft * (hxLeft + hxRight)) +
				beta * midDiff * (oldRightConc - oldLeftConc) *
					(rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight)) +
				beta * oldConc * (rightDiff - leftDiff) *
					(rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight));
		}
	}

	return;
}

bool
SoretDiffusionHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, double** concVector, double* val,
	int* indices, double hxLeft, double hxRight, int ix, double, int, double,
	int) const
{
	// Adjust the constants to the material
	double localHeatCond = getLocalHeatFactor(ix) * heatConductivity;

	if (ix == surfacePosition) {
		int diffClusterIdx = 0;

		for (auto const& currId : diffusingClusters) {
			auto cluster = network.getClusterCommon(currId);
			// Get the initial concentrations
			double oldConc = concVector[0][currId];
			double oldLeftConc = concVector[1][currId];
			double oldRightConc = concVector[2][currId];

			// Set the cluster index, the PetscSolver will use it to compute
			// the row and column indices for the Jacobian
			indices[diffClusterIdx] = currId;
			double midTemp = concVector[0][dof], rightTemp = concVector[2][dof];
			double midDiff = cluster.getDiffusionCoefficient(ix + 1),
				   rightDiff = cluster.getDiffusionCoefficient(ix + 2),
				   leftDiff = cluster.getDiffusionFactor() *
				exp(-cluster.getMigrationEnergy() /
					(kBoltzmann *
						(rightTemp +
							(hxLeft + hxRight) * heatFlux / localHeatCond)));
			// Compute the partial derivatives for diffusion of this cluster
			// for the middle, left, and right grid point
			val[diffClusterIdx * 6] =
				2.0 * beta * midDiff * heatFlux / (hxLeft * localHeatCond) -
				2.0 * midDiff / (hxLeft * hxRight) -
				beta * heatFlux * (rightDiff - leftDiff) /
					(localHeatCond * (hxLeft + hxRight)) -
				2.0 * beta * midDiff * heatFlux / (hxLeft * localHeatCond) -
				2.0 * beta * midDiff * (rightTemp - midTemp) /
					(hxLeft * hxRight) -
				beta * beta * midDiff * heatFlux * heatFlux /
					(localHeatCond * localHeatCond) +
				beta * heatFlux * (rightDiff - leftDiff) /
					(localHeatCond * (hxLeft + hxRight)); // middle conc
			val[(diffClusterIdx * 6) + 1] = 0.0; // left conc
			val[(diffClusterIdx * 6) + 2] =
				2.0 * midDiff / (hxLeft * hxRight); // right conc
			val[(diffClusterIdx * 6) + 3] = 2.0 * beta * midDiff * oldConc /
				(hxLeft * hxRight); // middle temp
			val[(diffClusterIdx * 6) + 4] = 0.0; // left temp
			val[(diffClusterIdx * 6) + 5] = -2.0 * beta * midDiff * oldConc /
				(hxLeft * hxRight); // right temp

			// Increase the index
			diffClusterIdx++;
		}
	}
	else {
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

			double leftTemp = concVector[1][dof], midTemp = concVector[0][dof],
				   rightTemp = concVector[2][dof];
			double leftDiff = cluster.getDiffusionCoefficient(ix),
				   midDiff = cluster.getDiffusionCoefficient(ix + 1),
				   rightDiff = cluster.getDiffusionCoefficient(ix + 2);

			// Compute the partial derivatives for diffusion of this cluster
			// for the middle, left, and right grid point
			val[diffClusterIdx * 6] = -2.0 * beta * midDiff *
					(leftTemp + (hxLeft / hxRight) * rightTemp -
						(1.0 + (hxLeft / hxRight)) * midTemp) /
					(hxLeft * (hxLeft + hxRight)) -
				beta * (rightDiff - leftDiff) * (rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight)); // middle conc
			val[(diffClusterIdx * 6) + 1] = beta * midDiff *
				(rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)); // left conc
			val[(diffClusterIdx * 6) + 2] = -beta * midDiff *
				(rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)); // right conc
			val[(diffClusterIdx * 6) + 3] = 2.0 * beta * midDiff * oldConc /
				(hxLeft * hxRight); // middle temp
			val[(diffClusterIdx * 6) + 4] = -2.0 * beta * midDiff * oldConc /
					(hxLeft * (hxLeft + hxRight)) +
				beta * midDiff * (oldRightConc - oldLeftConc) /
					((hxLeft + hxRight) * (hxLeft + hxRight)) +
				beta * oldConc * (rightDiff - leftDiff) /
					((hxLeft + hxRight) * (hxLeft + hxRight)); // left temp
			val[(diffClusterIdx * 6) + 5] = -2.0 * beta * midDiff * oldConc /
					(hxRight * (hxLeft + hxRight)) -
				beta * midDiff * (oldRightConc - oldLeftConc) /
					((hxLeft + hxRight) * (hxLeft + hxRight)) -
				beta * oldConc * (rightDiff - leftDiff) /
					((hxLeft + hxRight) * (hxLeft + hxRight)); // right temp

			// Increase the index
			diffClusterIdx++;
		}
	}

	return true;
}

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
