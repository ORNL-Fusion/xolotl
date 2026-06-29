#include <xolotl/core/network/FeReactionNetwork.h>
#include <xolotl/core/network/impl/FeReactionNetwork.tpp>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/util/Tokenizer.h>
#include <xolotl/core/network/impl/FeReaction.tpp>

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
template double //for SSBM connected to IReactionNetwork.h
ReactionNetwork<FeReactionNetwork>::getTotalVolumeRadius(
        ConcentrationsView concentrations, Species type, AmountType minSize);

template double //for SSBM connected to IReactionNetwork.h
ReactionNetwork<FeReactionNetwork>::getTotalRadiusVariance(
        ConcentrationsView concentrations, Species type, double mean ,AmountType minSize);

void //for SSBM 
FeReactionNetwork::initializeExtraClusterData(
	const options::IOptions& options)
{
	this->_clusterData.h_view().extraData.initialize(
		this->_clusterData.h_view().numClusters,
		this->_clusterData.h_view().gridSize);
	this->copyClusterDataView();

	auto data = this->_clusterData.h_view();
}

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

void
FeReactionNetwork::setReactionParams(std::string trapParams)
{
	// Convert the string to a vector of values
	auto tokens = util::Tokenizer<double>{trapParams}();
	// Set them in the corresponding reactions
	_reactions.forEachOn<FeTrapReaction>(
		"ReactionCollection::setReactionParams",
		DEVICE_LAMBDA(auto&& reaction) {
			auto i = reaction.getId();
			// Get the corresponding parameters
			reaction.setParameters(tokens[4 * i], tokens[4 * i + 1],
				tokens[4 * i + 2], tokens[4 * i + 3]);
		});
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
                           << speciesName << "_partial_density " << speciesName << "_partial_diameter ";
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
	
	double ahecomp_SSBMtotal = 0.0;
	double avcomp_SSBMtotal = 0.0;

	for (auto id = SpeciesId(numSpecies); id; ++id) {
		auto ms = minSizes[id()];
		auto totals = this->getTotals(conc,
			TQA{TQ{Q::total, id, 1}, TQ{Q::radius, id, 1}, TQ{Q::total, id, ms},
				TQ{Q::radius, id, ms}});

		totalVals[1 + (4 * id()) + 0] += totals[0] * fac;
		totalVals[1 + (4 * id()) + 1] += totals[1] * 2.0 * fac;
		totalVals[1 + (4 * id()) + 2] += totals[2] * fac;
		totalVals[1 + (4 * id()) + 3] += totals[3] * 2.0 * fac;
		
		// SSBM case
		if (this->_enableSSBM) {
			IndexType ssbmId = 0;
			IndexType ssbmSizeId = 0;

			// Compute average numbers
                        auto heConc = 0.0;
                        auto aheComp = 0.0;

                        auto vConc = 0.0;
                        auto avComp = 0.0;

                        switch (id()) {
                        // He
                        case 0:
                                ssbmId = this->_clusterData.h_view().bubbleId();
                                ssbmSizeId = ssbmId + 2;
                                vConc = conc(ssbmId);
                                if (vConc > 1.0e-16)
                                        avComp = conc(ssbmSizeId) / vConc;
					ahecomp_SSBMtotal += avComp * fac;
                                break;
                        // Void
                        case 1:
                                ssbmId = this->_clusterData.h_view().bubbleId();
                                ssbmSizeId = ssbmId + 1;
                                vConc = conc(ssbmId);
                                if (vConc > 1.0e-16)
                                        avComp = conc(ssbmSizeId) / vConc;
					avcomp_SSBMtotal += avComp * fac;
                                break;
                        default:
                                ssbmId = 0;
                                ssbmSizeId = 0;
                                break;
                        }

			// Add the single size data
			auto avRadius = 0.0;
			IndexType radiusId = 0;
			if (vConc > 1.0e-16) {
				if (id() > 4)
					radiusId = id() - 2;
				else if (id() > 1)
					radiusId = id() - 1;
				avRadius = util::max(0.0,
					computeBubbleRadius(avComp,
						this->_clusterData.h_view().latticeParameter()));
			}


			totalVals[1+(4 * id()) + 0] += vConc * fac;
			totalVals[1+(4 * id()) + 1] += vConc * avRadius * 2.0 * fac;
			
			if (avComp > minSizes[id()]) {
				totalVals[1+(4 * id()) + 2] += vConc * fac;
				totalVals[1+(4 * id()) + 3] += vConc * avRadius * 2.0 * fac;
			}


		}
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

			double totalHe = heConc * fac + ahecomp_SSBMtotal;
			double totalCav = cavConc +  avcomp_SSBMtotal;

			totalVals[0] += totalHe / totalCav;
			
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
			auto id = [i](std::size_t n) { return 1 + 4 * i + n; };
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
		outputFile << time << " " << globalData[0] << " ";
		for (auto i = 0; i < numSpecies; ++i) {
			auto id = [i](std::size_t n) { return 1 + 4 * i + n; };
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
