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
	double x = xGrid[xi + 1] - xGrid[surfacePosition + 1];
	if (x < interfaceLoc)
		return 0.2;
	return 1.0;
}

void
SoretDiffusionHandler::computeDiffusion(network::IReactionNetwork& network,
	double** concVector, double* updatedConcOffset, double hxLeft,
	double hxRight, int ix, double, int, double, int) const
{
	// Adjust the constants to the material
	double localHeatCond = getLocalHeatFactor(ix + localXs) * heatConductivity;

	// Surface
	if (ix + localXs == surfacePosition) {
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
	// Interface
	else if (fabs(xGrid[ix + localXs + 1] - xGrid[surfacePosition + 1] -
				 interfaceLoc) < 2.0) {
		double rightHeatCond =
			getLocalHeatFactor(ix + localXs + 1) * heatConductivity;
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
							(hxLeft + hxRight) * heatFlux / rightHeatCond)));

			updatedConcOffset[currId] += 2.0 *
					(J +
						(beta * midDiff * oldConc * heatFlux) / rightHeatCond) /
					hxLeft +
				2.0 * midDiff * (oldRightConc - oldConc) / (hxLeft * hxRight) -
				(J + (beta * midDiff * oldConc * heatFlux) / rightHeatCond) *
					(rightDiff - leftDiff) / (midDiff * (hxLeft + hxRight));

			// Second part
			updatedConcOffset[currId] -= 2.0 * beta * midDiff * oldConc *
					heatFlux / (hxLeft * rightHeatCond) +
				2.0 * beta * midDiff * oldConc * (rightTemp - midTemp) /
					(hxLeft * hxRight) +
				beta * heatFlux *
					(J +
						(beta * midDiff * oldConc * heatFlux) / rightHeatCond) /
					rightHeatCond -
				beta * oldConc * heatFlux * (rightDiff - leftDiff) /
					(rightHeatCond * (hxLeft + hxRight));
		}
	}
	// Bulk BC
	else if (fabs(xGrid[ix + localXs + 1] - xGrid[surfacePosition + 1] -
				 1000000) < 2.0) {
		// Neumann boundary for concentrations
		for (auto const& currId : diffusingClusters) {
			auto cluster = network.getClusterCommon(currId);
			// Get the initial concentrations
			double oldConc = concVector[0][currId];
			double oldLeftConc = concVector[1][currId];

			double leftTemp = concVector[1][dof], midTemp = concVector[0][dof],
				   rightTemp = concVector[2][dof];
			double leftDiff = cluster.getDiffusionCoefficient(ix),
				   midDiff = cluster.getDiffusionCoefficient(ix + 1),
				   rightDiff = cluster.getDiffusionCoefficient(ix + 2);

			updatedConcOffset[currId] +=
				2.0 * midDiff * (oldLeftConc - oldConc) / (hxLeft * hxRight) +
				2.0 * midDiff * beta * oldConc * (rightTemp - leftTemp) /
					(hxRight * (hxLeft + hxRight));

			// second part
			updatedConcOffset[currId] -= 2.0 * beta * midDiff * oldConc *
					(leftTemp + (hxLeft / hxRight) * rightTemp -
						(1.0 + (hxLeft / hxRight)) * midTemp) /
					(hxLeft * (hxLeft + hxRight)) +
				beta * midDiff * oldConc * beta * (rightTemp - leftTemp) *
					(rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight));
		}
	}
	else {
		for (auto const& currId : diffusingClusters) {
			auto cluster = network.getClusterCommon(currId);
			// Get the initial concentrations
			double oldConc = concVector[0][currId];
			double oldLeftConc = concVector[1][currId];
			double oldRightConc = concVector[2][currId];

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
	double localHeatCond = getLocalHeatFactor(ix + localXs) * heatConductivity;

	// Surface
	if (ix + localXs == surfacePosition) {
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
	// Interface
	else if (fabs(xGrid[ix + localXs + 1] - xGrid[surfacePosition + 1] -
				 interfaceLoc) < 2.0) {
		double rightHeatCond =
			getLocalHeatFactor(ix + localXs + 1) * heatConductivity;
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
							(hxLeft + hxRight) * heatFlux / rightHeatCond)));
			// Compute the partial derivatives for diffusion of this cluster
			// for the middle, left, and right grid point
			val[diffClusterIdx * 6] =
				2.0 * beta * midDiff * heatFlux / (hxLeft * rightHeatCond) -
				2.0 * midDiff / (hxLeft * hxRight) -
				beta * heatFlux * (rightDiff - leftDiff) /
					(rightHeatCond * (hxLeft + hxRight)) -
				2.0 * beta * midDiff * heatFlux / (hxLeft * rightHeatCond) -
				2.0 * beta * midDiff * (rightTemp - midTemp) /
					(hxLeft * hxRight) -
				beta * beta * midDiff * heatFlux * heatFlux /
					(rightHeatCond * rightHeatCond) +
				beta * heatFlux * (rightDiff - leftDiff) /
					(rightHeatCond * (hxLeft + hxRight)); // middle conc
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
	// Bulk BC
	else if (fabs(xGrid[ix + localXs + 1] - xGrid[surfacePosition + 1] -
				 1000000) < 2.0) {
		// Neumann boundary for concentrations
		int diffClusterIdx = 0;

		for (auto const& currId : diffusingClusters) {
			auto cluster = network.getClusterCommon(currId);
			// Set the cluster index, the PetscSolver will use it to compute
			// the row and column indices for the Jacobian
			indices[diffClusterIdx] = currId;
			// Get the initial concentrations
			double oldConc = concVector[0][currId];
			double oldLeftConc = concVector[1][currId];

			double leftTemp = concVector[1][dof], midTemp = concVector[0][dof],
				   rightTemp = concVector[2][dof];
			double leftDiff = cluster.getDiffusionCoefficient(ix),
				   midDiff = cluster.getDiffusionCoefficient(ix + 1),
				   rightDiff = cluster.getDiffusionCoefficient(ix + 2);

			// Compute the partial derivatives for diffusion of this cluster
			// for the middle, left, and right grid point
			val[diffClusterIdx * 6] = -2.0 * midDiff / (hxLeft * hxRight) +
				2.0 * beta * midDiff * (rightTemp - leftTemp) /
					(hxRight * (hxLeft + hxRight)) -
				2.0 * beta * midDiff *
					(leftTemp + (hxLeft / hxRight) * rightTemp -
						(1.0 + (hxLeft / hxRight)) * midTemp) /
					(hxLeft * (hxLeft + hxRight)) -
				beta * beta * midDiff * (rightTemp - leftTemp) *
					(rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight)); // middle conc
			val[(diffClusterIdx * 6) + 1] =
				2.0 * midDiff / (hxLeft * hxRight); // left conc
			val[(diffClusterIdx * 6) + 2] = 0.0; // right conc
			val[(diffClusterIdx * 6) + 3] = 2.0 * beta * midDiff * oldConc /
				(hxLeft * hxRight); // middle temp
			val[(diffClusterIdx * 6) + 4] = -2.0 * beta * midDiff * oldConc /
					(hxRight * (hxLeft + hxRight)) -
				2.0 * beta * midDiff * oldConc / (hxLeft * (hxLeft + hxRight)) -
				2.0 * beta * beta * midDiff * oldConc * (leftTemp - rightTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight)); // left temp
			val[(diffClusterIdx * 6) + 5] = -2.0 * beta * beta * midDiff *
				oldConc * (rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)); // right temp

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
