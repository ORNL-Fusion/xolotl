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
NEReactionNetwork::readReactions(double temperature, const std::string filename)
{
	syncClusterDataOnHost();

	// Read the reactions from a file
	std::ifstream reactionFile;
	reactionFile.open(filename);
	if (!reactionFile.good()) {
		this->_reactionEnergies =
			Kokkos::View<double***>("reactionEnergies", 0, 0, 0);
		return;
	}
	auto clData = _clusterData.h_view;

	auto nClusters = this->_numClusters;
	this->_reactionEnergies =
		Kokkos::View<double***>("reactionEnergies", nClusters, nClusters + 1,
			nClusters +
				1); // In case the second cluster or product has invalid index
	auto reactionEnergies = Kokkos::create_mirror_view(this->_reactionEnergies);
	auto constantRates =
		Kokkos::create_mirror_view(clData().extraData.constantRates);

	// The first loop will be on single clusters to get their properties
	std::vector<double> mVector(nClusters, 0.0);
	std::vector<double> g0Vector(nClusters, 0.0);
	std::vector<std::vector<std::pair<double, double>>> lVector(
		nClusters, std::vector<std::pair<double, double>>());
	constexpr double k_B = ::xolotl::core::kBoltzmann;

	// Get the line
	std::string line;
	getline(reactionFile, line);
	// Read the first line
	std::vector<double> tokens;
	util::Tokenizer<double>{line}(tokens);
	// And start looping on the lines
	while (tokens.size() > 0) {
		// Find the Id of the cluster
		Composition comp = Composition::zero();
		comp[Species::Xe] = static_cast<IndexType>(tokens[0]);
		comp[Species::V] = static_cast<IndexType>(tokens[1]);
		comp[Species::I] = static_cast<IndexType>(tokens[2]);
		auto rId = findCluster(comp, plsm::HostMemSpace{}).getId();

		// Get its properties and save it
		auto h0 = tokens[3];
		auto s0 = tokens[4];
		auto g0 = h0 - k_B * temperature * s0;
		auto D0 = tokens[5];
		auto q = tokens[6];
		auto mobility =
			D0 * exp(-q / (k_B * temperature)) / (k_B * temperature);
		g0Vector[rId] = g0;
		mVector[rId] = mobility;
		// Compute the thermal equilibrium
		auto theq = exp(-g0 / (k_B * temperature));

		// Loop on the linked clusters
		for (auto i = 7; i < tokens.size(); i += 4) {
			// Get its properties
			auto g0Linked = tokens[i] - k_B * temperature * tokens[i + 1];
			auto mobilityLinked = tokens[i + 2] *
				exp(-tokens[i + 3] / (k_B * temperature)) / (k_B * temperature);
			// Compute its thermal equilibrium
			auto theqLinked = exp(-g0Linked / (k_B * temperature));
			// Save the linkage information
			std::pair<double, double> linkage(
				mobilityLinked, theqLinked / theq);
			lVector[rId].push_back(linkage);
		}

		getline(reactionFile, line);
		if (line == "Reactions")
			break;

		tokens = util::Tokenizer<double>{line}();
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

		// Sink case
		if (tokens.size() == 4) {
			// Get the coefficient rate
			auto coefRate = tokens[3];
			// Compute the full given rate for this sink
			auto linkageRate = mVector[r1Id];
			// Loop on the linked clusters
			for (auto link : lVector[r1Id]) {
				linkageRate += link.first * link.second;
			}
			// Save the value
			constantRates(r1Id, nClusters) =
				coefRate * k_B * temperature * linkageRate;
			// Save the energy for this reaction
			reactionEnergies(r1Id, nClusters, nClusters) = -g0Vector[r1Id];
		}
		// Production case
		else {
			comp[Species::Xe] = static_cast<IndexType>(tokens[3]);
			comp[Species::V] = static_cast<IndexType>(tokens[4]);
			comp[Species::I] = static_cast<IndexType>(tokens[5]);
			auto r2Id = findCluster(comp, plsm::HostMemSpace{}).getId();

			comp[Species::Xe] = static_cast<IndexType>(tokens[6]);
			comp[Species::V] = static_cast<IndexType>(tokens[7]);
			comp[Species::I] = static_cast<IndexType>(tokens[8]);
			auto prodId = findCluster(comp, plsm::HostMemSpace{}).getId();

			// Get the coefficient rate
			auto coefRate = tokens[9];
			// Compute the base rate for this reaction
			auto linkageRate = mVector[r1Id] + mVector[r2Id];
			// Save the value to the left
			auto totalRate = coefRate * k_B * temperature * linkageRate;
			// Loop on the linked clusters
			for (auto link2 : lVector[r2Id]) {
				linkageRate += (link2.first + mVector[r1Id]) * link2.second;
				for (auto link1 : lVector[r1Id]) {
					linkageRate += (link2.first) * link2.second * link1.second;
				}
			}
			//			for (auto link2 : lVector[r2Id]) {
			//				linkageRate += (mVector[r1Id] + link2.first) *
			// link2.second;
			//			}
			// Save the value to the right
			totalRate = coefRate * k_B * temperature * linkageRate;
			constantRates(r1Id, r2Id) = totalRate;
			constantRates(r2Id, r1Id) = totalRate;

			// Save the energy for this reaction
			double energy = -g0Vector[r1Id] - g0Vector[r2Id];
			if (prodId != this->invalidIndex())
				energy += g0Vector[prodId];
			reactionEnergies(r1Id, r2Id, 0) = energy;
			reactionEnergies(r2Id, r1Id, 0) = energy;
		}

		getline(reactionFile, line);
		tokens.clear();
		util::Tokenizer<double>{line}(tokens);
	}

	deep_copy(this->_reactionEnergies, reactionEnergies);
	deep_copy(clData().extraData.constantRates, constantRates);
}

void
NEReactionNetwork::initializeExtraClusterData(const options::IOptions& options)
{
	if (!this->_enableReadRates) {
		return;
	}

	this->_clusterData.h_view().extraData.initialize(
		this->_clusterData.h_view().numClusters);
	this->copyClusterDataView();
}

namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
NEReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
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

	//	// Special case for Xe_1 + Xe_1V_1
	//	if (cl1Reg.isSimplex() && cl2Reg.isSimplex() &&
	//		((lo1.isOnAxis(Species::Xe) && lo1[Species::Xe] == 1 &&
	//			 lo2[Species::Xe] == 1 && lo2[Species::V] == 1) ||
	//			(lo2.isOnAxis(Species::Xe) && lo2[Species::Xe] == 1 &&
	//				lo1[Species::Xe] == 1 && lo1[Species::V] == 1))) {
	//		return;
	//	}

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
			this->addProductionReaction(tag, {i, j, k});
			// Loop on the species
			bool isOnAxis1 = false, isOnAxis2 = false;
			for (auto l : species) {
				if (lo1.isOnAxis(l()) && lo1[l()] == 1)
					isOnAxis1 = true;
				if (lo2.isOnAxis(l()) && lo2[l()] == 1)
					isOnAxis2 = true;
			}
			if (isOnAxis1 || isOnAxis2) {
				if (lo1.isOnAxis(Species::Xe) && lo2.isOnAxis(Species::Xe) &&
					lo1[Species::Xe] == 1 && lo2[Species::Xe] == 1) {
					continue;
				}

				this->addDissociationReaction(tag, {k, i, j});
			}
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
	if (clReg.isSimplex() && lo.isOnAxis(Species::V)) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}
}

inline ReactionCollection<NEReactionGenerator::NetworkType>
NEReactionGenerator::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->_clusterData.numClusters, this->_enableReadRates,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getReSolutionReactions(), this->getNucleationReactions(),
		this->getSinkReactions());
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
