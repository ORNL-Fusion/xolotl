#include <xolotl/core/network/ZrReactionNetwork.h>
#include <xolotl/core/network/impl/ZrReactionNetwork.tpp>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<ZrReactionNetwork>::ReactionNetwork();

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
ZrReactionNetwork::setConstantConnectivities(ConnectivitiesPair conns)
{
	_constantConnsRows = ConnectivitiesPairView(
		"dConstantConnectivitiesRows", conns.first.size());
	_constantConnsEntries = ConnectivitiesPairView(
		"dConstantConnectivitiesEntries", conns.second.size());
	auto hConnsRowsView = create_mirror_view(_constantConnsRows);
	auto hConnsEntriesView = create_mirror_view(_constantConnsEntries);
	for (auto i = 0; i < conns.first.size(); i++) {
		hConnsRowsView(i) = conns.first[i];
	}
	for (auto i = 0; i < conns.second.size(); i++) {
		hConnsEntriesView(i) = conns.second[i];
	}
	deep_copy(_constantConnsRows, hConnsRowsView);
	deep_copy(_constantConnsEntries, hConnsEntriesView);
}

void
ZrReactionNetwork::setConstantRateEntries()
{
	auto rows = _constantConnsRows;
	auto entries = _constantConnsEntries;
	_reactions.forEachOn<ZrConstantReaction>(
		"ReactionCollection::setConstantRates", DEVICE_LAMBDA(auto&& reaction) {
			reaction.defineRateEntries(rows, entries);
		});
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
					if (lo[Species::Basal] < data.transitionSize())
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

std::string
ZrReactionNetwork::getMonitorDataHeaderString() const
{
	std::stringstream header;

	auto numSpecies = getSpeciesListSize();
	header << "#time ";
	for (auto id = SpeciesId(numSpecies); id; ++id) {
		auto speciesName = this->getSpeciesName(id);
		header << speciesName << "_density " << speciesName << "_atom "
			   << speciesName << "_diameter " << speciesName
			   << "_partial_density " << speciesName << "_partial_atom "
			   << speciesName << "_partial_diameter ";
	}

	return header.str();
}

void
ZrReactionNetwork::addMonitorDataValues(Kokkos::View<const double*> conc,
	double fac, std::vector<double>& totalVals)
{
	auto numSpecies = getSpeciesListSize();
	const auto& minSizes = this->getMinRadiusSizes();
	for (auto id = SpeciesId(numSpecies); id; ++id) {
		using TQ = IReactionNetwork::TotalQuantity;
		using Q = TQ::Type;
		using TQA = util::Array<TQ, 6>;
		auto ms = minSizes[id()];
		auto totals = this->getTotals(conc,
			TQA{TQ{Q::total, id, 1}, TQ{Q::atom, id, 1}, TQ{Q::radius, id, 1},
				TQ{Q::total, id, ms}, TQ{Q::atom, id, ms},
				TQ{Q::radius, id, ms}});

		totalVals[(6 * id()) + 0] += totals[0] * fac;
		totalVals[(6 * id()) + 1] += totals[1] * fac;
		totalVals[(6 * id()) + 2] += totals[2] * 2.0 * fac;
		totalVals[(6 * id()) + 3] += totals[3] * fac;
		totalVals[(6 * id()) + 4] += totals[4] * fac;
		totalVals[(6 * id()) + 5] += totals[5] * 2.0 * fac;
	}
}

void
ZrReactionNetwork::writeMonitorDataLine(
	const std::vector<double>& localData, double time)
{
	auto numSpecies = getSpeciesListSize();

	// Sum all the concentrations through MPI reduce
	auto globalData = std::vector<double>(localData.size(), 0.0);
	MPI_Reduce(localData.data(), globalData.data(), localData.size(),
		MPI_DOUBLE, MPI_SUM, 0, util::getMPIComm());

	if (util::getMPIRank() == 0) {
		// Average the data
		for (auto i = 0; i < numSpecies; ++i) {
			auto id = [i](std::size_t n) { return 6 * i + n; };
			if (globalData[id(0)] > 1.0e-16) {
				globalData[id(2)] /= globalData[id(0)];
			}
			if (globalData[id(3)] > 1.0e-16) {
				globalData[id(5)] /= globalData[id(3)];
			}
		}

		// Set the output precision
		const int outputPrecision = 5;

		// Open the output file
		std::fstream outputFile;
		outputFile.open(
			getMonitorOutputFileName(), std::fstream::out | std::fstream::app);
		outputFile << std::setprecision(outputPrecision);

		// Output the data
		outputFile << time << " ";
		for (auto i = 0; i < numSpecies; ++i) {
			auto id = [i](std::size_t n) { return 6 * i + n; };
			outputFile << globalData[id(0)] << " " << globalData[id(1)] << " "
					   << globalData[id(2)] << " " << globalData[id(3)] << " "
					   << globalData[id(4)] << " " << globalData[id(5)] << " ";
		}
		outputFile << std::endl;

		// Close the output file
		outputFile.close();
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
