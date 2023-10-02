#include <xolotl/core/network/ZrReactionNetwork.h>
#include <xolotl/core/network/impl/ZrReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<ZrReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts);

template ReactionNetwork<ZrReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<ZrReactionNetwork>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);

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

			// adding basal
			auto size = hi[Species::V] + hi[Species::I] + hi[Species::Basal];

			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

void
ZrReactionNetwork::setConstantRates(RatesView rates, IndexType gridIndex)
{
	_reactions.forEachOn<ZrConstantReaction>(
		"ReactionCollection::setConstantRates", DEVICE_LAMBDA(auto&& reaction) {
			reaction.setRate(rates, gridIndex);
			reaction.updateRates();
		});
}

void
ZrReactionNetwork::setConstantConnectivities(ConnectivitiesVector conns)
{
	_constantConns = ConnectivitiesView(
		"dConstantConnectivities", conns.size(), conns[0].size());
	auto hConnsView = create_mirror_view(_constantConns);
	for (auto i = 0; i < conns.size(); i++)
		for (auto j = 0; j < conns[0].size(); j++) {
			hConnsView(i, j) = conns[i][j];
		}
	deep_copy(_constantConns, hConnsView);
}

void
ZrReactionNetwork::setGridSize(IndexType gridSize)
{
	this->_clusterData.h_view().extraData.setGridSize(
		this->_clusterData.h_view().numClusters, gridSize);
	Superclass::setGridSize(gridSize);
}

void
ZrReactionNetwork::initializeExtraClusterData(const options::IOptions& options)
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
			if (lo.isOnAxis(Species::V)) {
				// Spontaneous radii:
				// data.extraData.dislocationCaptureRadius(i, 0) = 3.05 *
				// pow(lo[Species::V], 0.12) / 10;
				// data.extraData.dislocationCaptureRadius(i, 1) = 0.39 *
				// pow(lo[Species::V], 0.4) / 10;

				// Thermal radii:
				if (lo[Species::V] < 1000) {
					data.extraData.dislocationCaptureRadius(i, 0) =
						2.8 * pow(lo[Species::V], 0.15) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						2.0 * pow(lo[Species::V], 0.3) / 10;
				}
				else {
					data.extraData.dislocationCaptureRadius(i, 0) = 0.79;
					data.extraData.dislocationCaptureRadius(i, 1) = 1.59;
				}
			}

			// adding basal
			// Set the dislocation capture radii for vacancy c-loops (convert to
			// nm): First index in dislocation capture radius is for I capture;
			// second is for V capture
			if (lo.isOnAxis(Species::Basal)) {
				// Spontaneous radii:
				// if(lo[Species::Basal] < ::xolotl::core::basalTransitionSize)
				// data.extraData.dislocationCaptureRadius(i, 0) = 3.9 *
				// pow(lo[Species::Basal], 0.07) / 10; if(lo[Species::Basal] <
				// ::xolotl::core::basalTransitionSize)
				// data.extraData.dislocationCaptureRadius(i, 1) = 0.55 *
				// pow(lo[Species::Basal], 0.33) / 10;

				// Thermal radii:
				if (lo[Species::Basal] < 1000) {
					if (lo[Species::Basal] <
						::xolotl::core::basalTransitionSize)
						data.extraData.dislocationCaptureRadius(i, 0) = 1.1;
					else
						data.extraData.dislocationCaptureRadius(i, 0) =
							5.2 * pow(lo[Species::Basal], 0.06) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						1.55 * pow(lo[Species::Basal], 0.28) / 10;
				}
				else {
					data.extraData.dislocationCaptureRadius(i, 0) = 0.787;
					data.extraData.dislocationCaptureRadius(i, 1) = 1.072;
				}
			}

			// Set the dislocation capture radii for interstitial a-loops
			// (convert to nm)
			else if (lo.isOnAxis(Species::I)) {
				// Spontaneous radii:
				// data.extraData.dislocationCaptureRadius(i, 0) = 4.2 *
				// pow(lo[Species::I], 0.05) / 10;
				// data.extraData.dislocationCaptureRadius(i, 1) = 5.1 *
				// pow(lo[Species::I], -0.01) / 10;

				// Thermal radii
				if (lo[Species::I] < 1000) {
					data.extraData.dislocationCaptureRadius(i, 0) =
						4.5 * pow(lo[Species::I], 0.205) / 10;
					data.extraData.dislocationCaptureRadius(i, 1) =
						6.0 * pow(lo[Species::I], 0.08) / 10;
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
