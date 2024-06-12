#include <xolotl/core/network/AlloyReactionNetwork.h>
#include <xolotl/core/network/impl/AlloyReactionNetwork.tpp>
#include <xolotl/util/MPIUtils.h>

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
			auto size =
				hi[Species::Void] + hi[Species::Frank] + hi[Species::Faulted];
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

void
AlloyReactionNetwork::setConstantRates(RatesView rates, IndexType gridIndex)
{
	_reactions.forEachOn<AlloyConstantReaction>(
		"ReactionCollection::setConstantRates", DEVICE_LAMBDA(auto&& reaction) {
			reaction.setRate(rates, gridIndex);
			reaction.updateRates();
		});
}

void
AlloyReactionNetwork::setConstantConnectivities(ConnectivitiesPair conns)
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
AlloyReactionNetwork::setConstantRateEntries()
{
	auto rows = _constantConnsRows;
	auto entries = _constantConnsEntries;
	_reactions.forEachOn<AlloyConstantReaction>(
		"ReactionCollection::setConstantRates", DEVICE_LAMBDA(auto&& reaction) {
			reaction.defineRateEntries(rows, entries);
		});
}

std::string
AlloyReactionNetwork::getMonitorDataHeaderString() const
{
	std::stringstream header;

	auto numSpecies = getSpeciesListSize();
	header << "#time ";
	for (auto id = SpeciesId(numSpecies); id; ++id) {
		auto speciesName = this->getSpeciesName(id);
		header << speciesName << "_density " << speciesName << "_diameter "
			   << speciesName << "_partial_density " << speciesName
			   << "_partial_diameter ";
	}

	return header.str();
}

void
AlloyReactionNetwork::addMonitorDataValues(Kokkos::View<const double*> conc,
	double fac, std::vector<double>& totalVals)
{
	auto numSpecies = getSpeciesListSize();
	const auto& minSizes = this->getMinRadiusSizes();
	for (auto id = SpeciesId(numSpecies); id; ++id) {
		using TQ = IReactionNetwork::TotalQuantity;
		using Q = TQ::Type;
		using TQA = util::Array<TQ, 4>;
		auto ms = minSizes[id()];
		auto totals = this->getTotals(conc,
			TQA{TQ{Q::total, id, 1}, TQ{Q::radius, id, 1}, TQ{Q::total, id, ms},
				TQ{Q::radius, id, ms}});

		totalVals[(4 * id()) + 0] += totals[0] * fac;
		totalVals[(4 * id()) + 1] += totals[1] * 2.0 * fac;
		totalVals[(4 * id()) + 2] += totals[2] * fac;
		totalVals[(4 * id()) + 3] += totals[3] * 2.0 * fac;
	}
}

void
AlloyReactionNetwork::writeMonitorDataLine(
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
			auto id = [i](std::size_t n) { return 4 * i + n; };
			if (globalData[id(0)] > 1.0e-16) {
				globalData[id(1)] /= globalData[id(0)];
			}
			if (globalData[id(2)] > 1.0e-16) {
				globalData[id(3)] /= globalData[id(2)];
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
			auto id = [i](std::size_t n) { return 4 * i + n; };
			outputFile << globalData[id(0)] << " " << globalData[id(1)] << " "
					   << globalData[id(2)] << " " << globalData[id(3)] << " ";
		}
		outputFile << std::endl;

		// Close the output file
		outputFile.close();
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
