#pragma once

#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>
#include <xolotl/core/network/impl/ZrClusterGenerator.tpp>
#include <xolotl/core/network/impl/ZrReaction.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
double
ZrReactionNetwork::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return alphaZrLatticeConstant;
	}
	return latticeParameter;
}

double
ZrReactionNetwork::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return alphaZrCoreRadius;
	}
	return impurityRadius;
}

ZrReactionNetwork::IndexType
ZrReactionNetwork::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = _clusterData.d_view;
	using Reducer = Kokkos::MaxLoc<ZrReactionNetwork::AmountType,
		ZrReactionNetwork::IndexType>;
	Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		_numClusters,
		KOKKOS_LAMBDA(IndexType i, Reducer::value_type & update) {
			const Region& clReg = clData().getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::V] + hi[Species::I];
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

void
ZrReactionNetwork::initializeExtraClusterData(const options::IOptions& options)
{
	if (!this->_enableSink) {
		return;
	}

	this->_clusterData.h_view().extraData.initialize(
		this->_clusterData.h_view().numClusters,
		this->_clusterData.h_view().gridSize);
	this->copyClusterDataView();

    auto data = this->_clusterData.h_view();
    Kokkos::parallel_for(
    this->_numClusters, KOKKOS_LAMBDA(const IndexType i) {
    auto cluster = data.getCluster(i);
    const auto& reg = cluster.getRegion();
    Composition lo(reg.getOrigin());

    // Set the dislocation capture radii for vacancy a-loops (convert to nm):
    // First index in dislocation capture radius is for I capture; second is for V capture
    if (lo.isOnAxis(Species::V)){
        data.extraData.dislocationCaptureRadius(i, 0) = 2.9 * pow(lo[Species::V], 0.12) / 10;
        data.extraData.dislocationCaptureRadius(i, 1) = 0.6 * pow(lo[Species::V], 0.3) / 10;
    }

    // Set the dislocation capture radii for interstitial a-loops (convert to nm)
    else if (lo.isOnAxis(Species::I)){
        data.extraData.dislocationCaptureRadius(i, 0) = 4.2 * pow(lo[Species::I], 0.05) / 10;
        data.extraData.dislocationCaptureRadius(i, 1) = 5.1 * pow(lo[Species::I], -0.01) / 10;
    }

    // ADD BASAL DISLOCATION CAPTURE RADII HERE

    });
}

namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
ZrReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	using Species = typename Network::Species;
	using Composition = typename Network::Composition;
	using AmountType = typename Network::AmountType;

	if (i == j) {
		addSinks(i, tag);
	}

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	// vac + vac = vac
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::V] + lo2[Species::V];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::V] = size;
		auto vProdId = subpaving.findTileId(comp, plsm::onDevice);
		if (vProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, vProdId});
			if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
				this->addDissociationReaction(tag, {vProdId, i, j});
			}
		}
		return;
	}

	// vac + int = vac | int | recombine
	if (((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
			(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I)))) {
		// Find out which one is which
		auto vSize =
			lo1.isOnAxis(Species::V) ? lo1[Species::V] : lo2[Species::V];
		auto iSize =
			lo1.isOnAxis(Species::I) ? lo1[Species::I] : lo2[Species::I];
		// Compute the product size
		int prodSize = vSize - iSize;
		// 3 cases
		if (prodSize > 0) {
			// Looking for V cluster
			Composition comp = Composition::zero();
			comp[Species::V] = prodSize;
			auto vProdId = subpaving.findTileId(comp, plsm::onDevice);
			if (vProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, vProdId});
				// No dissociation
			}
		}
		else if (prodSize < 0) {
			// Looking for I cluster
			Composition comp = Composition::zero();
			comp[Species::I] = -prodSize;
			auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
			if (iProdId != subpaving.invalidIndex()) {
				this->addProductionReaction(tag, {i, j, iProdId});
				// No dissociation
			}
		}
		else {
			// No product
			this->addProductionReaction(tag, {i, j});
		}
		return;
	}

	// int + int = int
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I)) {
		// Compute the composition of the new cluster
		auto size = lo1[Species::I] + lo2[Species::I];
		// Find the corresponding cluster
		Composition comp = Composition::zero();
		comp[Species::I] = size;
		auto iProdId = subpaving.findTileId(comp, plsm::onDevice);
		if (iProdId != subpaving.invalidIndex()) {
			this->addProductionReaction(tag, {i, j, iProdId});
			if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
				this->addDissociationReaction(tag, {iProdId, i, j});
			}
		}
		return;
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
ZrReactionGenerator::addSinks(IndexType i, TTag tag) const
{
	using Species = typename Network::Species;
	using Composition = typename Network::Composition;

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	// I
	if (clReg.isSimplex() && lo.isOnAxis(Species::I)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}

	// V
	if (clReg.isSimplex() && lo.isOnAxis(Species::V)) {
		this->addSinkReaction(tag, {i, Network::invalidIndex()});
	}
}

inline ReactionCollection<ZrReactionGenerator::Network>
ZrReactionGenerator::getReactionCollection() const
{
	ReactionCollection<Network> ret(this->_clusterData.gridSize,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions());
	return ret;
}
} // namespace detail

inline detail::ZrReactionGenerator
ZrReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::ZrReactionGenerator{*this};
}

namespace detail
{
KOKKOS_INLINE_FUNCTION
void
ZrClusterUpdater::updateDiffusionCoefficient(
	const ClusterData& data, IndexType clusterId, IndexType gridIndex) const
{
	// I migration energies in eV
	constexpr Kokkos::Array<double, 6> iMigrationA = {
		0.0, 0.17, 0.23, 0.49, 0.75, 0.87};
	constexpr Kokkos::Array<double, 6> iMigrationC = {
		0.0, 0.30, 0.54, 0.93, 1.2, 1.6};
    //I diffusion factors in nm^2/s
	constexpr Kokkos::Array<double, 6> iDiffusionA = {
		0.0, 2.4e+11, 3.2e+11, 4.9e+12, 5.1e+13, 4.3e+13};
	constexpr Kokkos::Array<double, 6> iDiffusionC = {
		0.0, 6.8e+11, 2.6e+12, 6.8e+13, 4.2e+14, 5.5e+15};


	// V migration energies in eV
	constexpr Kokkos::Array<double, 10> vMigrationA = {
		0.0, 0.59, 0.58, 0.94, 0.16, 0.81, 0.25, 0.90, 0.90, 0.90};
    constexpr Kokkos::Array<double, 10> vMigrationC = {
		0.0, 0.67, 0.41, 1.12, 0.58, 0.29, 0.18, 1.0, 1.0, 1.0};
    // V diffusions factors in nm^2/s
	constexpr Kokkos::Array<double, 10> vDiffusionA = {
		0.0,1.6e+12, 2.7e+12, 4.9e+13, 2.5e+10, 2e+13, 3.2e+10, 1.5e+13, 1.5e+13, 1.5e+13};
	constexpr Kokkos::Array<double, 10> vDiffusionC = {
		0.0, 2.2e+12, 2.3e+11, 1.27e+15, 4.5e+11, 5.7e+11, 9.1e+9, 4.2e+13, 4.2e+13, 4.2e+13};

	// 3D diffuser case
	if (data.migrationEnergy(clusterId) < 0.0) {
		double kernel = -1.0 / (kBoltzmann * data.temperature(gridIndex));
		const auto& clReg = data.getCluster(clusterId).getRegion();
		Network::Composition lo = clReg.getOrigin();
		using Species = Network::Species;

		if (lo.isOnAxis(Species::I)) {
			// Compute each contribution
			double Da = iDiffusionA[lo[Species::I]] *
				exp(iMigrationA[lo[Species::I]] * kernel);
			double Dc = iDiffusionC[lo[Species::I]] *
				exp(iMigrationC[lo[Species::I]] * kernel);

			// Compute the mean
			data.diffusionCoefficient(clusterId, gridIndex) =
				pow(Da * Da * Dc, 1.0 / 3.0);

			// Compute the anisotropy factor
			data.extraData.anisotropyRatio(clusterId, gridIndex) =
				pow(Dc / Da, 1.0 / 6.0);
			return;
		}

		if (lo.isOnAxis(Species::V)) {
			// Compute each contribution
			double Da = vDiffusionA[lo[Species::V]] *
				exp(vMigrationA[lo[Species::V]] * kernel);
			double Dc = vDiffusionC[lo[Species::V]] *
				exp(vMigrationC[lo[Species::V]] * kernel);

			// Compute the mean
			data.diffusionCoefficient(clusterId, gridIndex) =
				pow(Da * Da * Dc, 1.0 / 3.0);

			// Compute the anisotropy factor
			data.extraData.anisotropyRatio(clusterId, gridIndex) =
				pow(Dc / Da, 1.0 / 6.0);

			return;
		}
	}

	// 1D diffuser case
	data.diffusionCoefficient(clusterId, gridIndex) =
		data.diffusionFactor(clusterId) *
		exp(-data.migrationEnergy(clusterId) /
			(kBoltzmann * data.temperature(gridIndex)));
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
