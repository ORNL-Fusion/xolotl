// Includes
#include <array>
#include <iostream>

#include <xolotl/core/Constants.h>
#include <xolotl/core/modified/SoretDiffusionHandler.h>
#include <xolotl/core/network/IPSIReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace modified
{
template <typename T>
using HostUnmanaged =
	Kokkos::View<T, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

void
SoretDiffusionHandler::syncDiffusingClusters(network::IReactionNetwork& network)
{
	auto clusterIds_h = HostUnmanaged<const IdType*>(
		diffusingClusters.data(), diffusingClusters.size());
	diffClusterIds = Kokkos::View<IdType*>(
		"Diffusing Cluster Ids", diffusingClusters.size());
	deep_copy(diffClusterIds, clusterIds_h);

	diffClusters = Kokkos::View<DeviceCluster*>(
		Kokkos::ViewAllocateWithoutInitializing("Diffusing Clusters"),
		diffusingClusters.size());
	auto clusters_h = create_mirror_view(diffClusters);
	for (IdType i = 0; i < diffusingClusters.size(); ++i) {
		clusters_h[i] = network.getClusterCommon(
			diffusingClusters[i], plsm::DeviceMemSpace{});
	}
	deep_copy(diffClusters, clusters_h);

	auto beta_h =
		HostUnmanaged<const double*>(betaFactor.data(), betaFactor.size());
	beta = Kokkos::View<double*>(
		Kokkos::ViewAllocateWithoutInitializing("beta factor"),
		betaFactor.size());
	deep_copy(beta, beta_h);
}

void
SoretDiffusionHandler::initialize(network::IReactionNetwork& network,
	std::vector<core::RowColPair>& idPairs, std::vector<double> grid, int xs)
{
	// Clear the index vector
	diffusingClusters.clear();

	// Set dof
	dof = network.getDOF();

	// Copy the grid
	xGrid = grid;
	localXs = xs;

	using NetworkType = network::IPSIReactionNetwork;
	using AmountType = NetworkType::AmountType;

	auto psiNetwork = dynamic_cast<NetworkType*>(&network);
	if (!psiNetwork) {
		throw std::runtime_error(
			"Soret Diffusion only supported for PSI reaction network types");
	}
	auto numSpecies = psiNetwork->getSpeciesListSize();
	auto specIdHe = psiNetwork->getHeliumSpeciesId();

	// Initialize the composition
	auto comp = std::vector<AmountType>(numSpecies, 0);

	auto addCluster = [&](const std::vector<AmountType>& comp,
						  double betaValue) {
		// Check that the cluster is present in the network
		auto clusterId = psiNetwork->findClusterId(comp);
		if (clusterId != NetworkType::invalidIndex()) {
			auto cluster = network.getClusterCommon(clusterId);

			// Get its diffusion factor and migration energy
			double diffFactor = cluster.getDiffusionFactor();

			// Don't do anything if the diffusion factor is 0.0
			if (not util::equal(diffFactor, 0.0)) {
				// Note that cluster is diffusing.
				diffusingClusters.emplace_back(clusterId);
				betaFactor.emplace_back(betaValue);

				// This cluster interacts with temperature now
				idPairs.push_back(core::RowColPair{clusterId, dof});
			}
		}
	};

	// Helium
	comp[specIdHe()] = 1;
	addCluster(comp, 0.0065);
	comp[specIdHe()] = 0;

	// Hydrogen
	if (psiNetwork->hasDeuterium()) {
		auto clusterSpecies = network.parseSpeciesId("D");
		comp[clusterSpecies()] = 1;
		addCluster(comp, 0.0045);
		comp[clusterSpecies()] = 0;
	}
	if (psiNetwork->hasTritium()) {
		auto clusterSpecies = network.parseSpeciesId("T");
		comp[clusterSpecies()] = 1;
		addCluster(comp, 0.0045);
		comp[clusterSpecies()] = 0;
	}

	// Self-interstitial
	auto specIdI = psiNetwork->getInterstitialSpeciesId();
	comp[specIdI()] = 1;
	addCluster(comp, 0.0128);

	this->syncDiffusingClusters(network);
}

void
SoretDiffusionHandler::computeDiffusion(network::IReactionNetwork& network,
	const StencilConcArray& concVector, Kokkos::View<double*> updatedConcOffset,
	double hxLeft, double hxRight, int ix, double, int, double, int) const
{
	if (concVector.size() != 3) {
		throw std::runtime_error(
			"Wrong size for 1D concentration stencil; should be 3, got " +
			std::to_string(concVector.size()));
	}

	Kokkos::Array<Kokkos::View<const double*>, 3> concVec = {
		concVector[0], concVector[1], concVector[2]};

	auto clusterIds = this->diffClusterIds;
	auto clusters = this->diffClusters;
	auto beta = this->beta;
	Kokkos::parallel_for(
		clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
			auto currId = clusterIds[i];
			auto cluster = clusters[i];

			// Get the initial concentrations
			double oldConc = concVec[0][currId];
			double oldLeftConc = concVec[1][currId];
			double oldRightConc = concVec[2][currId];
			double leftDiff = cluster.getDiffusionCoefficient(ix);
			double midDiff = cluster.getDiffusionCoefficient(ix + 1);
			double rightDiff = cluster.getDiffusionCoefficient(ix + 2);
			double leftTemp = cluster.getTemperature(ix);
			double midTemp = cluster.getTemperature(ix + 1);
			double rightTemp = cluster.getTemperature(ix + 2);

			// Update the concentration of the cluster
			// TODO: Should this use atomic_?
			updatedConcOffset[currId] -= 2.0 * beta[i] * midDiff * oldConc *
					(leftTemp + (hxLeft / hxRight) * rightTemp -
						(1.0 + (hxLeft / hxRight)) * midTemp) /
					(hxLeft * (hxLeft + hxRight)) +
				beta[i] * midDiff * (oldRightConc - oldLeftConc) *
					(rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight)) +
				beta[i] * oldConc * (rightDiff - leftDiff) *
					(rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight));
		});

	// TODO: Maybe we need a Kokkos::fence() here?
}

bool
SoretDiffusionHandler::computePartialsForDiffusion(
	network::IReactionNetwork& network, const StencilConcArray& concVector,
	Kokkos::View<double*> values, double hxLeft, double hxRight, int ix, double,
	int, double, int) const
{
	Kokkos::Array<Kokkos::View<const double*>, 3> concVec = {
		concVector[0], concVector[1], concVector[2]};

	auto clusterIds = this->diffClusterIds;
	auto clusters = this->diffClusters;
	auto beta = this->beta;
	Kokkos::parallel_for(
		clusterIds.size(), KOKKOS_LAMBDA(IdType i) {
			auto currId = clusterIds[i];
			auto cluster = clusters[i];

			// Get the initial concentrations
			double oldConc = concVec[0][currId];
			double oldLeftConc = concVec[1][currId];
			double oldRightConc = concVec[2][currId];
			double leftDiff = cluster.getDiffusionCoefficient(ix);
			double midDiff = cluster.getDiffusionCoefficient(ix + 1);
			double rightDiff = cluster.getDiffusionCoefficient(ix + 2);
			double leftTemp = cluster.getTemperature(ix);
			double midTemp = cluster.getTemperature(ix + 1);
			double rightTemp = cluster.getTemperature(ix + 2);

			// Compute the partial derivatives for diffusion of this cluster
			// for the middle, left, and right grid point
			// TODO: Should this use atomic_?
			values[i * 3] = -2.0 * beta[i] * midDiff *
					(leftTemp + (hxLeft / hxRight) * rightTemp -
						(1.0 + (hxLeft / hxRight)) * midTemp) /
					(hxLeft * (hxLeft + hxRight)) -
				beta[i] * (rightDiff - leftDiff) * (rightTemp - leftTemp) /
					((hxLeft + hxRight) * (hxLeft + hxRight)); // middle conc

			values[(i * 3) + 1] = beta[i] * midDiff * (rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)); // left conc

			values[(i * 3) + 2] = -beta[i] * midDiff * (rightTemp - leftTemp) /
				((hxLeft + hxRight) * (hxLeft + hxRight)); // right conc
		});

	// TODO: Maybe we need a Kokkos::fence() here?

	return true;
}

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
