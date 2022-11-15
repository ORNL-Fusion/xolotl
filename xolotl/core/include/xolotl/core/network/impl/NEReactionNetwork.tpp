#pragma once

#include <xolotl/core/network/detail/impl/NucleationReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/ReSolutionReactionGenerator.tpp>
#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/impl/NEClusterGenerator.tpp>
#include <xolotl/core/network/impl/NEReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
void
NEReactionNetwork::readClusters(const std::string filename)
{
	// Read the reactions from a file
	std::ifstream reactionFile;
	reactionFile.open(filename);
	if (!reactionFile.good()) {
		return;
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
	fileClusterNumber = 0;
	while (tokens.size() > 0) {
		fileClusterNumber++;

		getline(reactionFile, line);
		if (line == "Reactions")
			break;

		tokens = util::Tokenizer<double>{line}();
	}
}

void
NEReactionNetwork::readReactions(double temperature, const std::string filename)
{
	syncClusterDataOnHost();

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
	std::vector<double> mVector(fileClusterNumber, 0.0);
	std::vector<double> g0Vector(fileClusterNumber, 0.0);
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
		auto mobility =
			D0 * exp(-q / (k_B * temperature)) / (k_B * temperature);
		g0Vector[count] = g0;
		mVector[count] = mobility;
		// Save the linkage information
		std::pair<double, double> linkage(mobility, g0);
		lVector[count].push_back(linkage);

		// Loop on the linked clusters
		for (auto i = 7; i < tokens.size(); i += 4) {
			// Get its properties
			auto g0Linked = tokens[i] - k_B * temperature * tokens[i + 1];
			auto mobilityLinked = tokens[i + 2] *
				exp(-tokens[i + 3] / (k_B * temperature)) / (k_B * temperature);
			// Save the linkage information
			linkage.first = mobilityLinked;
			linkage.second = g0Linked;
			lVector[count].push_back(linkage);
		}
		// Set values in the cluster
		if (_clusterDataMirror.has_value()) {
			// Compute the diffusivity with linkage
			double diffusivity = 0.0;
			double largestG0 = g0Vector[count];
			IdType largestId = 0;
			for (auto link : lVector[count]) {
				diffusivity += link.first *
					exp((largestG0 - link.second) / (k_B * temperature));
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

			double largestG0 = g0Vector[map1Id];
			IdType largestId = 0;

			// Compute the full given rate for this sink
			auto linkageRate = lVector[map1Id][largestId].first;
			constantRates(map1Id, fileClusterNumber, 1) =
				coefRate * k_B * temperature * linkageRate;
			// Loop on the linked clusters
			linkageRate = 0.0;
			for (auto link : lVector[map1Id]) {
				linkageRate += link.first *
					exp((largestG0 - link.second) / (k_B * temperature));
			}
			// Save the value
			constantRates(map1Id, fileClusterNumber, 0) =
				coefRate * k_B * temperature * linkageRate;

			double energy = 0.0;
			for (auto i = 0; i < 1; i++) {
				auto g0 = lVector[map1Id][i].second;
				double li = exp((largestG0 - g0) / (k_B * temperature));
				energy -= li * g0;
			}
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
			auto mapProdId =
				fileClusterMap.value_at(fileClusterMap.find(prodId));

			double largestG01 = g0Vector[map1Id];
			IdType largestId1 = 0;
			double largestG02 = g0Vector[map2Id];
			IdType largestId2 = 0;

			// Get the coefficient rate
			//			double coefRate = static_cast<double>(tokens[9]);
			double coefRate = 1.0;
			// Compute the base rate for this reaction
			auto linkageRate = 0.0;
			// Loop on the linked clusters
			for (auto i = 0; i < lVector[map1Id].size(); i++) {
				for (auto j = 0; j < lVector[map2Id].size(); j++) {
					double d_a = lVector[map1Id][i].first;
					double l_a = exp((largestG01 - lVector[map1Id][i].second) /
						(k_B * temperature));
					double d_b = lVector[map2Id][j].first;
					double l_b = exp((largestG02 - lVector[map2Id][j].second) /
						(k_B * temperature));
					linkageRate += (d_a + d_b) * l_a * l_b;
				}
			}
			// Save the value to the right
			auto totalRate = coefRate * k_B * temperature * linkageRate;
			constantRates(map1Id, map2Id, 0) = totalRate;
			constantRates(map2Id, map1Id, 0) = totalRate;

			// Save the energy for this reaction
			double energy = 0.0;
			for (auto i = 0; i < 1; i++) {
				auto g0 = lVector[map1Id][i].second;
				double li = exp((largestG01 - g0) / (k_B * temperature));
				energy -= li * g0;
			}
			for (auto i = 0; i < 1; i++) {
				auto g0 = lVector[map2Id][i].second;
				double li = exp((largestG02 - g0) / (k_B * temperature));
				energy -= li * g0;
			}

			if (prodId != this->invalidIndex()) {
				double largestG0 = g0Vector[mapProdId];
				IdType largestId = 0;
				//				for (auto i = 0; i < lVector[prodId].size();
				// i++)
				//{ 					if (lVector[prodId][i].second < largestG0)
				//{ largestG0 =
				// lVector[prodId][i].second; 						largestId =
				// i;
				//					}
				//				}

				linkageRate = 0.0;
				for (auto i = 0; i < lVector[map1Id].size(); i++) {
					for (auto j = 0; j < lVector[map2Id].size(); j++) {
						for (auto k = 0; k < lVector[mapProdId].size(); k++) {
							double d_a = lVector[map1Id][i].first;
							double d_b = lVector[map2Id][j].first;
							double l_z =
								exp((largestG0 - lVector[mapProdId][k].second) /
									(k_B * temperature));
							linkageRate += (d_a + d_b) * l_z;
						}
					}
				}

				totalRate = coefRate * k_B * temperature * linkageRate;
				constantRates(map1Id, map2Id, 1) = totalRate;
				constantRates(map2Id, map1Id, 1) = totalRate;

				// Add linkage
				for (auto i = 0; i < 1; i++) {
					auto g0 = lVector[mapProdId][i].second;
					double li = exp((largestG0 - g0) / (k_B * temperature));
					energy += li * g0;
				}
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

	invalidateDataMirror();
	copyClusterDataView();
}

void
NEReactionNetwork::initializeExtraClusterData(const options::IOptions& options)
{
	if (!this->_enableReadRates) {
		return;
	}

	this->_clusterData.h_view().extraData.initialize(fileClusterNumber);
	this->copyClusterDataView();
}

namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
NEReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	// Check the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;
	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}

	// TODO: re-solution and nucleation need to be added
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();
	constexpr auto invalidIndex = NetworkType::invalidIndex();

	auto numClusters = this->getNumberOfClusters();

	if (i == j) {
		addSinks(i, tag);
	}

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	auto& subpaving = this->getSubpaving();

	// Special case for I + V
	if (cl1Reg.isSimplex() && cl2Reg.isSimplex() &&
		((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
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
			auto vProdId = subpaving.findTileId(comp);
			if (vProdId != invalidIndex) {
				this->addProductionReaction(tag, {i, j, vProdId});
				this->addDissociationReaction(tag, {vProdId, i, j});
			}
		}
		else if (prodSize < 0) {
			// Looking for I cluster
			Composition comp = Composition::zero();
			comp[Species::I] = -prodSize;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != invalidIndex) {
				this->addProductionReaction(tag, {i, j, iProdId});
				this->addDissociationReaction(tag, {iProdId, i, j});
			}
		}
		else {
			// No product
			this->addProductionReaction(tag, {i, j});
		}
		return;
	}

	// General case
	constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
	using BoundsArray =
		Kokkos::Array<Kokkos::pair<AmountType, AmountType>, numSpeciesNoI>;
	plsm::EnumIndexed<BoundsArray, Species> bounds;
	// Loop on the species
	for (auto l : species) {
		auto low = lo1[l] + lo2[l];
		auto high = hi1[l] + hi2[l] - 2;
		// Special case for I
		if (l == Species::I) {
			bounds[Species::V].first -= high;
			bounds[Species::V].second -= low;
		}
		else {
			bounds[l] = {low, high};
		}
	}

	// Look for potential product
	IndexType nProd = 0;
	for (IndexType k = 0; k < numClusters; ++k) {
		// Get the composition
		const auto& prodReg = this->getCluster(k).getRegion();
		bool isGood = true;
		// Loop on the species
		for (auto l : speciesNoI) {
			if (prodReg[l()].begin() > bounds[l()].second) {
				isGood = false;
				break;
			}
			if (prodReg[l()].end() - 1 < bounds[l()].first) {
				isGood = false;
				break;
			}
		}

		if (isGood) {
			nProd++;
			this->addProductionReaction(tag, {i, j, k});
			this->addDissociationReaction(tag, {k, i, j});
		}
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
NEReactionGenerator::addSinks(IndexType i, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	const auto& clReg = this->getCluster(i).getRegion();
	Composition lo = clReg.getOrigin();

	// I
	if (clReg.isSimplex() && lo.isOnAxis(Species::I)) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}

	// V
	if (clReg.isSimplex() && lo.isOnAxis(Species::V) && lo[Species::V] < 3) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}
}

inline ReactionCollection<NEReactionGenerator::NetworkType>
NEReactionGenerator::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->_clusterData.extraData.fileClusterMap.size(),
		this->_enableReadRates, this->getProductionReactions(),
		this->getDissociationReactions(), this->getReSolutionReactions(),
		this->getNucleationReactions(), this->getSinkReactions());
	ret._data.reactionEnergies = this->_reactionEnergies;
	return ret;
}
} // namespace detail

inline detail::NEReactionGenerator
NEReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::NEReactionGenerator{*this};
}

namespace detail
{
KOKKOS_INLINE_FUNCTION
void
NEClusterUpdater::updateDiffusionCoefficient(
	const ClusterData& data, IndexType clusterId, IndexType gridIndex) const
{
	if (data.migrationEnergy(clusterId) > 0.0) {
		//		using Species = typename NetworkType::Species;
		//		using Composition = typename NetworkType::Composition;
		//		const auto& clReg = data.tiles(clusterId).getRegion();
		//		Composition lo = clReg.getOrigin();
		//		if (clReg.isSimplex() && lo[Species::V] == 1 && lo[Species::Xe]
		//== 1) {
		//			// Intrinsic diffusion
		//			double kernel = -3.04 / (kBoltzmann *
		// data.temperature(gridIndex)); 			double D3 = 7.6e8 *
		// exp(kernel);
		// // nm2/s
		//
		//			// We need the fission rate now
		//			double fissionRate = data.fissionRate(0) * 1.0e27; // #/m3/s
		//
		//			// Athermal diffusion
		//			double D1 = (8e-40 * fissionRate) * 1.0e18; // nm2/s
		//
		//			// Radiation-enhanced diffusion
		//			kernel = -1.2 / (kBoltzmann * data.temperature(gridIndex));
		//			double D2 =
		//				(5.6e-25 * sqrt(fissionRate) * exp(kernel)) * 1.0e18; //
		// nm2/s
		//
		//			data.diffusionCoefficient(clusterId, gridIndex) = D1 + D2 +
		// D3;
		//
		//			return;
		//		}

		data.diffusionCoefficient(clusterId, gridIndex) =
			data.diffusionFactor(clusterId) *
			exp(-data.migrationEnergy(clusterId) /
				(kBoltzmann * data.temperature(gridIndex)));

		return;
	}

	// If the diffusivity is given
	data.diffusionCoefficient(clusterId, gridIndex) =
		data.diffusionFactor(clusterId);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
