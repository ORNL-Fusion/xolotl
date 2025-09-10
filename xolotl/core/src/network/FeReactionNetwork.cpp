#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/impl/FeReactionNetwork.tpp>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<FeReactionNetwork>::ReactionNetwork();

template ReactionNetwork<FeReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts);

template ReactionNetwork<FeReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts);

template double
ReactionNetwork<FeReactionNetwork>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeReactionNetwork>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeReactionNetwork>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeReactionNetwork>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<FeReactionNetwork>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);

double
FeReactionNetwork::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return ironLatticeConstant;
	}
	return latticeParameter;
}

double
FeReactionNetwork::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return heliumRadius;
	}
	return impurityRadius;
}

FeReactionNetwork::IndexType
FeReactionNetwork::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = _clusterData.d_view;
	using Reducer = Kokkos::MaxLoc<FeReactionNetwork::AmountType,
		FeReactionNetwork::IndexType>;
	Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		"FeReactionNetwork::checkLargestClusterId", _numClusters,
		KOKKOS_LAMBDA(IndexType i, Reducer::value_type & update) {
			const Region& clReg = clData().getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::He] + hi[Species::V];
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

std::string
FeReactionNetwork::getMonitorDataHeaderString() const
{
	std::stringstream header;

	auto numSpecies = getSpeciesListSize();
	header << "#time He_cav ";
	for (auto id = SpeciesId(numSpecies); id; ++id) {
		auto speciesName = this->getSpeciesName(id);
		header << speciesName << "_density " << speciesName << "_diameter "
			   << speciesName << "_partial_density " << speciesName
			   << "_partial_diameter ";
	}

	return header.str();
}

void
FeReactionNetwork::addMonitorDataValues(Kokkos::View<const double*> conc,
	double fac, std::vector<double>& totalVals)
{
	auto numSpecies = getSpeciesListSize();
	const auto& minSizes = this->getMinRadiusSizes();
	using TQ = IReactionNetwork::TotalQuantity;
	using Q = TQ::Type;
	using TQA = util::Array<TQ, 4>;
	for (auto id = SpeciesId(numSpecies); id; ++id) {
		auto ms = minSizes[id()];
		auto totals = this->getTotals(conc,
			TQA{TQ{Q::total, id, 1}, TQ{Q::radius, id, 1}, TQ{Q::total, id, ms},
				TQ{Q::radius, id, ms}});

		totalVals[1 + (4 * id()) + 0] += totals[0] * fac;
		totalVals[1 + (4 * id()) + 1] += totals[1] * 2.0 * fac;
		totalVals[1 + (4 * id()) + 2] += totals[2] * fac;
		totalVals[1 + (4 * id()) + 3] += totals[3] * 2.0 * fac;

		// Special case for trapped helium
		if (id() == 0) {
			// Find the vacancy index
			constexpr auto speciesRangeNoI = getSpeciesRangeNoI();
			bool hasVacancy = false;
			Species vIndex;
			for (auto i : speciesRangeNoI) {
				if (isVacancy(i)) {
					hasVacancy = true;
					vIndex = i;
				}
			}

			auto tiles = _subpaving.getTiles();
			double heConc = 0.0;
			Kokkos::parallel_reduce(
				"FeReactionNetwork::TrappedAtom", this->_numClusters,
				KOKKOS_LAMBDA(IndexType i, double& lsum) {
					const Region& clReg = tiles(i).getRegion();
					if (clReg[vIndex].begin() > 0) {
						const auto factor = clReg.volume() / clReg[id].length();
						for (AmountType j : makeIntervalRange(clReg[id])) {
							if (j >= 1)
								lsum += conc(i) * j * factor;
						}
					}
				},
				heConc);
			double cavConc = 0.0;
			Kokkos::parallel_reduce(
				"FeReactionNetwork::TrappedAtom", this->_numClusters,
				KOKKOS_LAMBDA(IndexType i, double& lsum) {
					const Region& clReg = tiles(i).getRegion();
					if (clReg[vIndex].begin() > 0) {
						const auto factor = clReg.volume() / clReg[id].length();
						for (AmountType j : makeIntervalRange(clReg[id])) {
							if (j >= 1)
								lsum += conc(i) * factor;
						}
					}
				},
				cavConc);

			Kokkos::fence();

			totalVals[0] += heConc * fac / cavConc;
		}
	}
}

void
FeReactionNetwork::writeMonitorDataLine(
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
