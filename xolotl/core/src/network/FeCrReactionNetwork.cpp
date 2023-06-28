#include <xolotl/core/network/FeCrReactionNetwork.h>
#include <xolotl/core/network/impl/FeCrReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<FeCrReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts);

template ReactionNetwork<FeCrReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts);

template double
ReactionNetwork<FeCrReactionNetwork>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeCrReactionNetwork>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeCrReactionNetwork>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeCrReactionNetwork>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeCrReactionNetwork>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeCrReactionNetwork>::getSmallConcentration(
	ConcentrationsView concentrations, Species type, AmountType maxSize);

double
FeCrReactionNetwork::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return ironLatticeConstant;
	}
	return latticeParameter;
}

double
FeCrReactionNetwork::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return fecrCoreRadius;
	}
	return impurityRadius;
}

FeCrReactionNetwork::IndexType
FeCrReactionNetwork::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = _clusterData.d_view;
	using Reducer = Kokkos::MaxLoc<FeCrReactionNetwork::AmountType,
		FeCrReactionNetwork::IndexType>;
	Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		"FeCrReactionNetwork::checkLargestClusterId", _numClusters,
		KOKKOS_LAMBDA(IndexType i, Reducer::value_type & update) {
			const Region& clReg = clData().getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::Free] + hi[Species::Trapped] +
				hi[Species::Junction] + hi[Species::Loop];
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}
  
void
FeCrReactionNetwork::initializeExtraClusterData(
	const options::IOptions& options)
{
	this->_clusterData.h_view().extraData.initialize(
		this->_clusterData.h_view().numClusters);
	this->copyClusterDataView();
	this->invalidateDataMirror();
}

void
FeCrReactionNetwork::computeFluxesPreProcess(ConcentrationsView concentrations,
	FluxesView fluxes, IndexType gridIndex, double surfaceDepth, double spacing)
{
	auto data = this->_clusterData.h_view();
	auto diffusionFactor = this->getClusterDataMirror().diffusionFactor;
	auto sigma = create_mirror_view(data.extraData.netSigma);
	for (auto i = 0; i < this->_numClusters; i++) {
		//		sigma(i) = 0.02;
		if (diffusionFactor(i) == 0.0) {
			sigma(i) = 0.0;
			continue;
		}
		sigma(i) = this->getNetSigma(concentrations, i, gridIndex);
	}
	deep_copy(data.extraData.netSigma, sigma);
	this->updateReactionRates(0.0);
}

void
FeCrReactionNetwork::computePartialsPreProcess(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex, double surfaceDepth, double spacing)
{
	auto data = this->_clusterData.h_view();
	auto diffusionFactor = this->getClusterDataMirror().diffusionFactor;
	auto sigma = create_mirror_view(data.extraData.netSigma);
	for (auto i = 0; i < this->_numClusters; i++) {
		//		sigma(i) = 0.02;
		if (diffusionFactor(i) == 0.0) {
			sigma(i) = 0.0;
			continue;
		}
		sigma(i) = this->getNetSigma(concentrations, i, gridIndex);
	}
	deep_copy(data.extraData.netSigma, sigma);
	this->updateReactionRates(0.0);
}

void
FeCrReactionNetwork::updateOutgoingSinkFluxes(
	double* gridPointSolution, std::vector<double>& fluxes, IndexType gridIndex)
{
	auto clusterData = this->getClusterDataMirror();
	auto density = this->getClusterDataMirror().sinkDensity(); // nm-2
	auto portion =
		this->getClusterDataMirror().sinkPortion(); // Portion of screw
	auto r = 1.0 / sqrt(::xolotl::core::pi * density); // nm
	auto rCore = ::xolotl::core::fecrCoreRadius;
	auto temperature = this->getClusterDataMirror().temperature(0);
	constexpr double K = 170.0e9; // GPa
	constexpr double nu = 0.29;
	constexpr double b = 0.25; // nm
	double deltaV =
		1.67 * this->getClusterDataMirror().atomicVolume() * 1.0e-27; // m3
	//	constexpr double a0 = 0.91, a1 = -2.16, a2 = -0.92; // Random dipole
	constexpr double a0 = 0.87, a1 = -5.12, a2 = -0.77; // Full network
	constexpr double k_B = 1.380649e-23; // J K-1.

	double L = (K * b * deltaV * (1.0 - 2.0 * nu)) /
		(2.0 * ::xolotl::core::pi * k_B * temperature * (1.0 - nu));

	double delta = sqrt(rCore * rCore + (L * L) / 4.0);

	double edge = (2.0 * ::xolotl::core::pi * (a0 + a1 * (rCore / r)) *

					  (1.0 - portion) *
					  ((std::log(r / rCore) *
						   (a0 * r + a1 * delta + a2 * (delta - rCore))) /
						  (std::log(r / delta) * (a0 * r + a1 * rCore)))) /
		(std::log(r / rCore));

	double screw =
		(2.0 * ::xolotl::core::pi * (a0 + a1 * (rCore / r)) * portion) /
		(std::log(r / rCore));

	auto diffusionFactor = this->getClusterDataMirror().diffusionFactor;
	auto dof = this->getDOF();
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(gridPointSolution, dof);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof);
	deep_copy(dConcs, hConcs);

	// Loop on all the clusters
	for (auto i = 0; i < this->_numClusters; ++i) {
		if (diffusionFactor(i) == 0.0)
			continue;

		auto cluster = this->getCluster(i, plsm::HostMemSpace{});
		const auto& clReg = cluster.getRegion();
		Composition lo = clReg.getOrigin();
		auto diffCoef = cluster.getDiffusionCoefficient(gridIndex);
		auto factor = density * gridPointSolution[i] * diffCoef;

		// V case
		if (lo[Species::V] > 0) {
			double size =
				lo[Species::V] + (double)(clReg[Species::V].length() - 1) / 2.0;
			// edge
			fluxes[0] += factor * edge * size;
			// screw
			fluxes[1] += factor * screw * size;
		}
		else {
			// I and Free case
			double bias = 1.0;
			double size = lo[Species::Free] +
				(double)(clReg[Species::Free].length() - 1) / 2.0;
			if (lo[Species::I] > 0) {
				bias = 1.05;
				size = lo[Species::I] +
					(double)(clReg[Species::I].length() - 1) / 2.0;
				double r0 = (r + rCore);
				edge *= bias;
			}
			else {
				auto sigma = this->getNetSigma(dConcs, i, gridIndex);
				edge *= sigma;
				screw *= sigma;
			}
			// edge
			fluxes[2] += factor * edge * size;
			// screw
			fluxes[3] += factor * screw * size;
		}
	}

	return;
}
} // namespace network
} // namespace core
} // namespace xolotl
