#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/impl/AlloyReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<AlloyReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts);

template ReactionNetwork<AlloyReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts);

template double
ReactionNetwork<AlloyReactionNetwork>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<AlloyReactionNetwork>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<AlloyReactionNetwork>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<AlloyReactionNetwork>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<AlloyReactionNetwork>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);

double
AlloyReactionNetwork::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return alloyLatticeConstant;
	}
	return latticeParameter;
}

double
AlloyReactionNetwork::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return alloyCoreRadius;
	}
	return impurityRadius;
}

AlloyReactionNetwork::IndexType
AlloyReactionNetwork::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = _clusterData.d_view;
	using Reducer = Kokkos::MaxLoc<AlloyReactionNetwork::AmountType,
		AlloyReactionNetwork::IndexType>;
	Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		"AlloyReactionNetwork::checkLargestClusterId", _numClusters,
		KOKKOS_LAMBDA(IndexType i, Reducer::value_type & update) {
			const Region& clReg = clData().getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::PerfectV] + hi[Species::FaultedI] +
				hi[Species::FaultedV] + hi[Species::PerfectI];
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

void
AlloyReactionNetwork::initializeExtraClusterData(
	const options::IOptions& options)
{
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

			// Set the dislocation capture radii for vacancy a-loops (convert to
			// nm): First index in dislocation capture radius is for I capture;
			// second is for V capture
			if (lo.isOnAxis(Species::PerfectV)) {
				// Thermal radii:
				if (lo[Species::PerfectV] < 1000) {
					data.extraData.dislocationCaptureRadius(i, 0) =
						2.8 * pow(lo[Species::PerfectV], 0.15) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						2.0 * pow(lo[Species::PerfectV], 0.3) / 10;
				}
				else {
					data.extraData.dislocationCaptureRadius(i, 0) = 0.79;
					data.extraData.dislocationCaptureRadius(i, 1) = 1.59;
				}
			}
			else if (lo.isOnAxis(Species::FaultedV)) {
				// Thermal radii:
				if (lo[Species::FaultedV] < 1000) {
					data.extraData.dislocationCaptureRadius(i, 0) =
						2.8 * pow(lo[Species::FaultedV], 0.15) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						2.0 * pow(lo[Species::FaultedV], 0.3) / 10;
				}
				else {
					data.extraData.dislocationCaptureRadius(i, 0) = 0.79;
					data.extraData.dislocationCaptureRadius(i, 1) = 1.59;
				}
			}

			// Set the dislocation capture radii for interstitial a-loops
			// (convert to nm)
			else if (lo.isOnAxis(Species::PerfectI)) {
				// Thermal radii
				if (lo[Species::PerfectI] < 1000) {
					data.extraData.dislocationCaptureRadius(i, 0) =
						4.5 * pow(lo[Species::PerfectI], 0.205) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						6.0 * pow(lo[Species::PerfectI], 0.08) / 10;
				}
				else {
					data.extraData.dislocationCaptureRadius(i, 0) = 1.85;
					data.extraData.dislocationCaptureRadius(i, 1) = 1.04;
				}
			}

			// Set the dislocation capture radii for interstitial a-loops
			// (convert to nm)
			else if (lo.isOnAxis(Species::FaultedI)) {
				// Thermal radii
				if (lo[Species::FaultedI] < 1000) {
					data.extraData.dislocationCaptureRadius(i, 0) =
						4.5 * pow(lo[Species::FaultedI], 0.205) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						6.0 * pow(lo[Species::FaultedI], 0.08) / 10;
				}
				else {
					data.extraData.dislocationCaptureRadius(i, 0) = 1.85;
					data.extraData.dislocationCaptureRadius(i, 1) = 1.04;
				}
			}
		}); // Goes with parallel_for
}
} // namespace network
} // namespace core
} // namespace xolotl
