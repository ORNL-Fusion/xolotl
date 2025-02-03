#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/core/network/impl/NEReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
template ReactionNetwork<NEReactionNetwork>::ReactionNetwork();

template ReactionNetwork<NEReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts,
	const std::vector<SubdivisionRatio>& subdivisionRatios, IndexType gridSize,
	const options::IOptions& opts);

template ReactionNetwork<NEReactionNetwork>::ReactionNetwork(
	const std::vector<AmountType>& maxSpeciesAmounts, IndexType gridSize,
	const options::IOptions& opts);

template double
ReactionNetwork<NEReactionNetwork>::getTotalConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<NEReactionNetwork>::getTotalRadiusConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<NEReactionNetwork>::getTotalAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<NEReactionNetwork>::getTotalTrappedAtomConcentration(
	ConcentrationsView concentrations, Species type, AmountType minSize);

template double
ReactionNetwork<NEReactionNetwork>::getTotalVolumeFraction(
	ConcentrationsView concentrations, Species type, AmountType minSize);

double
NEReactionNetwork::checkLatticeParameter(double latticeParameter)
{
	if (latticeParameter <= 0.0) {
		return uraniumDioxydeLatticeConstant;
	}
	return latticeParameter;
}

double
NEReactionNetwork::checkImpurityRadius(double impurityRadius)
{
	if (impurityRadius <= 0.0) {
		return xenonRadius;
	}
	return impurityRadius;
}

NEReactionNetwork::IndexType
NEReactionNetwork::checkLargestClusterId()
{
	// Copy the cluster data for the parallel loop
	auto clData = _clusterData.d_view;
	using Reducer = Kokkos::MaxLoc<NEReactionNetwork::AmountType,
		NEReactionNetwork::IndexType>;
	Reducer::value_type maxLoc;
	Kokkos::parallel_reduce(
		"NEReactionNetwork::checkLargestClusterId", _numClusters,
		KOKKOS_LAMBDA(IndexType i, Reducer::value_type & update) {
			const Region& clReg = clData().getCluster(i).getRegion();
			Composition hi = clReg.getUpperLimitPoint();
			auto size = hi[Species::Xe] + hi[Species::V];
			if (size > update.val) {
				update.val = size;
				update.loc = i;
			}
		},
		Reducer(maxLoc));

	return maxLoc.loc;
}

void
NEReactionNetwork::initializeExtraClusterData(const options::IOptions& options)
{
	if (!this->_enableReadRates) {
		return;
	}

	this->_clusterData.h_view().extraData.initialize(
		this->_clusterData.h_view().extraData.fileClusterNumber);
	this->copyClusterDataView();
}

void
NEReactionNetwork::readClusters(const std::string filename)
{
	// Read the reactions from a file
	std::ifstream reactionFile;
	reactionFile.open(filename);
	if (!reactionFile.good()) {
		throw std::string("\nIncorrect filename for reading reactions! \n"
						  "Either provide reaction.dat or set the file path "
						  "with the reactionFilePath option.");
	}

	// First find how big of an object we need to create
	// i.e. how many clusters in the file
	// Get the line
	std::string line;
	getline(reactionFile, line);
	// Read the first line
	std::vector<double> tokens;
	util::Tokenizer<double>{line}(tokens);
	// And start looping on the lines
	IdType count = 0;
	while (tokens.size() > 0) {
		count++;

		getline(reactionFile, line);
		if (line == "Reactions")
			break;

		tokens = util::Tokenizer<double>{line}();
	}
	this->_clusterData.h_view().extraData.fileClusterNumber = count;
}

void
NEReactionNetwork::readReactions(double temperature, const std::string filename)
{
	syncClusterDataOnHost();

	auto fileClusterNumber =
		this->_clusterData.h_view().extraData.fileClusterNumber;

	// Read the reactions from a file
	std::ifstream reactionFile;
	reactionFile.open(filename);
	if (!reactionFile.good()) {
		this->_reactionEnergies =
			Kokkos::View<double**>("reactionEnergies", 0, 0);
		return;
	}
	auto clData = _clusterData.h_view;

	this->_reactionEnergies =
		Kokkos::View<double**>("reactionEnergies", fileClusterNumber,
			fileClusterNumber +
				1); // In case the second cluster or product has invalid index
	auto reactionEnergies = Kokkos::create_mirror_view(this->_reactionEnergies);
	auto constantRates =
		Kokkos::create_mirror_view(clData().extraData.constantRates);
	auto fileClusterMap =
		Kokkos::UnorderedMap<int, int, Kokkos::DefaultHostExecutionSpace>(
			fileClusterNumber);

	// The first loop will be on single clusters to get their properties
	std::vector<double> g0Vector(fileClusterNumber, 0.0);
	std::vector<IdType> smallestG0Id(fileClusterNumber, 0);
	std::vector<std::vector<std::pair<double, double>>> lVector(
		fileClusterNumber, std::vector<std::pair<double, double>>());
	constexpr double k_B = ::xolotl::core::kBoltzmann;

	// Get the line
	std::string line;
	getline(reactionFile, line);
	// Read the first line
	std::vector<double> tokens;
	util::Tokenizer<double>{line}(tokens);
	// And start looping on the lines
	int count = 0;
	while (tokens.size() > 0) {
		// Find the Id of the cluster
		Composition comp = Composition::zero();
		comp[Species::Xe] = static_cast<IndexType>(tokens[0]);
		comp[Species::V] = static_cast<IndexType>(tokens[1]);
		comp[Species::I] = static_cast<IndexType>(tokens[2]);

		auto rId = findCluster(comp, plsm::HostMemSpace{}).getId();

		fileClusterMap.insert(rId, count);

		// Get its properties and save it
		auto h0 = tokens[3];
		auto s0 = tokens[4];
		auto g0 = h0 - k_B * temperature * s0;
		auto D0 = tokens[5];
		auto q = tokens[6];
		auto diffusivity = D0 * exp(-q / (k_B * temperature));
		// Save the linkage information
		std::pair<double, double> linkage(diffusivity, g0);
		lVector[count].push_back(linkage);
		g0Vector[count] = g0;

		// Loop on the linked clusters
		for (auto i = 7; i < tokens.size(); i += 4) {
			// Get its properties
			auto g0Linked = tokens[i] - k_B * temperature * tokens[i + 1];
			auto diffLinked =
				tokens[i + 2] * exp(-tokens[i + 3] / (k_B * temperature));
			// Save the linkage information
			linkage.first = diffLinked;
			linkage.second = g0Linked;
			lVector[count].push_back(linkage);

			// Look for the smallest G0
			if (g0Linked < g0Vector[count]) {
				smallestG0Id[count] = (i - 7) / 4 + 1;
				g0Vector[count] = g0Linked;
			}
		}
		// Set values in the cluster
		if (_clusterDataMirror.has_value()) {
			// Compute the diffusivity with linkage
			double diffusivity = 0.0;
			for (auto link : lVector[count]) {
				diffusivity += link.first *
					exp((g0Vector[count] - link.second) / (k_B * temperature));
			}
			this->_clusterDataMirror.value().diffusionFactor(rId) = diffusivity;
			this->_clusterDataMirror.value().migrationEnergy(rId) = 0.0;
		}

		getline(reactionFile, line);
		if (line == "Reactions")
			break;

		tokens = util::Tokenizer<double>{line}();
		count++;
	}

	// Now we need to loop on the reactions to set their rates
	getline(reactionFile, line);
	// Read the next line
	tokens.clear();
	util::Tokenizer<double>{line}(tokens);
	// And start looping on the lines
	while (tokens.size() > 0) {
		// Find the Id of the first reactant
		Composition comp = Composition::zero();
		comp[Species::Xe] = static_cast<IndexType>(tokens[0]);
		comp[Species::V] = static_cast<IndexType>(tokens[1]);
		comp[Species::I] = static_cast<IndexType>(tokens[2]);
		auto r1Id = findCluster(comp, plsm::HostMemSpace{}).getId();
		auto map1Id = fileClusterMap.value_at(fileClusterMap.find(r1Id));

		// Sink case
		if (tokens.size() == 4) {
			// Get the coefficient rate
			auto coefRate = tokens[3];

			// Compute the full given rate for this sink
			auto linkageRate = lVector[map1Id][smallestG0Id[map1Id]].first;
			constantRates(map1Id, fileClusterNumber, 1) =
				coefRate * linkageRate;
			// Loop on the linked clusters
			linkageRate = 0.0;
			for (auto link : lVector[map1Id]) {
				linkageRate += link.first *
					exp((g0Vector[map1Id] - link.second) / (k_B * temperature));
			}
			// Save the value
			constantRates(map1Id, fileClusterNumber, 0) =
				coefRate * linkageRate;

			double energy = -g0Vector[map1Id];
			reactionEnergies(map1Id, fileClusterNumber) = energy;
		}
		// Production case
		else {
			comp[Species::Xe] = static_cast<IndexType>(tokens[3]);
			comp[Species::V] = static_cast<IndexType>(tokens[4]);
			comp[Species::I] = static_cast<IndexType>(tokens[5]);
			auto r2Id = findCluster(comp, plsm::HostMemSpace{}).getId();
			auto map2Id = fileClusterMap.value_at(fileClusterMap.find(r2Id));

			comp[Species::Xe] = static_cast<IndexType>(tokens[6]);
			comp[Species::V] = static_cast<IndexType>(tokens[7]);
			comp[Species::I] = static_cast<IndexType>(tokens[8]);
			auto prodId = findCluster(comp, plsm::HostMemSpace{}).getId();

			if (static_cast<int>(tokens[9]) == 0) {
				constantRates(map1Id, map2Id, 0) = 0.0;
				constantRates(map2Id, map1Id, 0) = 0.0;
				constantRates(map1Id, map2Id, 1) = 0.0;
				constantRates(map2Id, map1Id, 1) = 0.0;
				reactionEnergies(map1Id, map2Id) = 0.0;
				reactionEnergies(map2Id, map1Id) = 0.0;

				getline(reactionFile, line);
				tokens.clear();
				util::Tokenizer<double>{line}(tokens);

				continue;
			}

			// Get the coefficient rate
			//			double coefRate = static_cast<double>(tokens[9]);
			double coefRate = 1.0;
			// Compute the base rate for this reaction
			auto linkageRate = 0.0;
			// Loop on the linked clusters
			for (auto i = 0; i < lVector[map1Id].size(); i++) {
				for (auto j = 0; j < lVector[map2Id].size(); j++) {
					double d_a = lVector[map1Id][i].first;
					double l_a =
						exp((g0Vector[map1Id] - lVector[map1Id][i].second) /
							(k_B * temperature));
					double d_b = lVector[map2Id][j].first;
					double l_b =
						exp((g0Vector[map2Id] - lVector[map2Id][j].second) /
							(k_B * temperature));
					linkageRate += (d_a + d_b) * l_a * l_b;
				}
			}
			// Save the value to the right
			auto totalRate = coefRate * linkageRate;
			constantRates(map1Id, map2Id, 0) = totalRate;
			constantRates(map2Id, map1Id, 0) = totalRate;

			// Save the energy for this reaction
			double energy = -g0Vector[map1Id] - g0Vector[map2Id];

			if (prodId != this->invalidIndex()) {
				auto mapProdId =
					fileClusterMap.value_at(fileClusterMap.find(prodId));
				linkageRate = 0.0;
				for (auto i = 0; i < lVector[map1Id].size(); i++) {
					for (auto j = 0; j < lVector[map2Id].size(); j++) {
						for (auto k = 0; k < lVector[mapProdId].size(); k++) {
							double d_a = lVector[map1Id][i].first;
							double d_b = lVector[map2Id][j].first;
							double l_z = exp((g0Vector[mapProdId] -
												 lVector[mapProdId][k].second) /
								(k_B * temperature));
							linkageRate += (d_a + d_b) * l_z;
						}
					}
				}

				totalRate = coefRate * linkageRate;
				constantRates(map1Id, map2Id, 1) = totalRate;
				constantRates(map2Id, map1Id, 1) = totalRate;

				// Add energy
				energy += g0Vector[mapProdId];
			}
			reactionEnergies(map1Id, map2Id) = energy;
			reactionEnergies(map2Id, map1Id) = energy;
		}

		getline(reactionFile, line);
		tokens.clear();
		util::Tokenizer<double>{line}(tokens);
	}

	// Find the Id of the interstitial cluster
	Composition comp = Composition::zero();
	comp[Species::I] = 1;
	auto rId = findCluster(comp, plsm::HostMemSpace{}).getId();
	auto mapId = fileClusterMap.value_at(fileClusterMap.find(rId));
	// Save its formation energy
	clData().setIFormationEnergy(g0Vector[mapId]);

	// Same with vacancy
	comp[Species::I] = 0;
	comp[Species::V] = 1;
	rId = findCluster(comp, plsm::HostMemSpace{}).getId();
	mapId = fileClusterMap.value_at(fileClusterMap.find(rId));
	// Save its formation energy
	clData().setVFormationEnergy(g0Vector[mapId]);
	comp[Species::V] = 2;
	rId = findCluster(comp, plsm::HostMemSpace{}).getId();
	mapId = fileClusterMap.value_at(fileClusterMap.find(rId));
	// Save its formation energy
	clData().setV2FormationEnergy(g0Vector[mapId]);

	// Same with xenon
	comp[Species::V] = 0;
	comp[Species::Xe] = 1;
	rId = findCluster(comp, plsm::HostMemSpace{}).getId();
	mapId = fileClusterMap.value_at(fileClusterMap.find(rId));
	// Save its formation energy
	clData().setXeFormationEnergy(g0Vector[mapId]);

	deep_copy(this->_reactionEnergies, reactionEnergies);
	deep_copy(clData().extraData.constantRates, constantRates);
	deep_copy(clData().extraData.fileClusterMap, fileClusterMap);
	deep_copy(clData().diffusionFactor,
		this->_clusterDataMirror.value().diffusionFactor);
	deep_copy(clData().migrationEnergy,
		this->_clusterDataMirror.value().migrationEnergy);

	this->invalidateDataMirror();
	this->copyClusterDataView();
}
} // namespace network
} // namespace core
} // namespace xolotl
