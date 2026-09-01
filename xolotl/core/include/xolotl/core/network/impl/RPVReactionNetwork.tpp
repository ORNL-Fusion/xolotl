#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/impl/SinkReactionGenerator.tpp>
#include <xolotl/core/network/impl/RPVClusterGenerator.tpp>
#include <xolotl/core/network/impl/RPVReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
void
RPVReactionNetwork::computeFluxesPreProcess(ConcentrationsView concentrations,
	FluxesView fluxes, IndexType gridIndex, double surfaceDepth, double spacing)
{
	// Set S diffusivity from V concentration
	// Get the single vacancy cluster
	auto& subpaving = this->getSubpaving();
	Composition comp = Composition::zero();
	comp[Species::V] = 1;
	auto vId = subpaving.findTileId(comp);
	if (vId == subpaving.invalidIndex())
		throw std::runtime_error(
			"The Vacancy cluster is not present in the network");

	auto singleVConc = concentrations(vId);

	// Compute the equilibrium concentration with a 1.85 eV formation energy
	double kernel =
		-1.0 / (kBoltzmann * _clusterData.h_view().temperature(gridIndex));
	auto concEq = exp(kernel * 1.85) / _atomicVolume;

	// Get the solute cluster
	comp[Species::V] = 0;
	comp[Species::S] = 1;
	auto sId = subpaving.findTileId(comp);
	if (sId == subpaving.invalidIndex())
		throw std::runtime_error(
			"The Solute cluster is not present in the network");

	// Get its migration energy and diffusion factor
	auto migration = _clusterData.h_view().migrationEnergy(sId);
	auto diffFactor = _clusterData.h_view().diffusionFactor(sId);

	// Set the diffusion coefficient
	double diffCoef =
		diffFactor * exp(kernel * migration) * singleVConc / concEq;
	// double diffCoef = 0.0;
	_clusterData.h_view().diffusionCoefficient(sId, gridIndex) = diffCoef;

	// Set sink term from the precipitate population
	// Get the radius concentration R * C above a given size
	auto radiusConc =
		getTotalRadiusConcentration(concentrations, Species::S, 44);
	_clusterData.h_view().setPrecipitateStrength(
		4.0 * ::xolotl::core::pi * radiusConc);

	// Update the sink rates
	using SinkReactionType = typename Superclass::Traits::SinkReactionType;
	auto sinkReactions = this->_reactions.template getView<SinkReactionType>();
	Kokkos::parallel_for(
		"RPVReactionNetwork::updateReactionRates", sinkReactions.size(),
		KOKKOS_LAMBDA(IndexType i) { sinkReactions[i].updateRates(); });

	invalidateDataMirror();
}

void
RPVReactionNetwork::computePartialsPreProcess(ConcentrationsView concentrations,
	Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
	double spacing)
{
	// Set S diffusivity from V concentration
	// Get the single vacancy cluster
	auto& subpaving = this->getSubpaving();
	Composition comp = Composition::zero();
	comp[Species::V] = 1;
	auto vId = subpaving.findTileId(comp);
	if (vId == subpaving.invalidIndex())
		throw std::runtime_error(
			"The Vacancy cluster is not present in the network");

	auto singleVConc = concentrations(vId);

	// Compute the equilibrium concentration with a 1.85 eV formation energy
	double kernel =
		-1.0 / (kBoltzmann * _clusterData.h_view().temperature(gridIndex));
	auto concEq = exp(kernel * 1.85) / _atomicVolume;

	// Get the solute cluster
	comp[Species::V] = 0;
	comp[Species::S] = 1;
	auto sId = subpaving.findTileId(comp);
	if (sId == subpaving.invalidIndex())
		throw std::runtime_error(
			"The Solute cluster is not present in the network");

	// Get its migration energy and diffusion factor
	auto migration = _clusterData.h_view().migrationEnergy(sId);
	auto diffFactor = _clusterData.h_view().diffusionFactor(sId);

	// Set the diffusion coefficient
	double diffCoef =
		diffFactor * exp(kernel * migration) * singleVConc / concEq;
	// double diffCoef = 0.0;
	_clusterData.h_view().diffusionCoefficient(sId, gridIndex) = diffCoef;

	// Set sink term from the precipitate population
	// Get the radius concentration R * C above a given size
	auto radiusConc =
		getTotalRadiusConcentration(concentrations, Species::S, 44);
	_clusterData.h_view().setPrecipitateStrength(
		4.0 * ::xolotl::core::pi * radiusConc);

	// Update the sink rates
	using SinkReactionType = typename Superclass::Traits::SinkReactionType;
	auto sinkReactions = this->_reactions.template getView<SinkReactionType>();
	Kokkos::parallel_for(
		"RPVReactionNetwork::updateReactionRates", sinkReactions.size(),
		KOKKOS_LAMBDA(IndexType i) { sinkReactions[i].updateRates(); });

	invalidateDataMirror();
}

namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
RPVReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	// Check the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;
	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}

	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	if (i == j) {
		addSinks(i, tag);
	}

	auto numClusters = this->getNumberOfClusters();

	auto& subpaving = this->getSubpaving();
	auto previousIndex = subpaving.invalidIndex();

	// I + I = I and Loop
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::I)) {
		// Compute the composition of the new cluster
		auto minSize = lo1[Species::I] + lo2[Species::I];
		auto maxSize = hi1[Species::I] + hi2[Species::I] - 2;

		Composition comp = Composition::zero();
		// Loop on the possibilities
		for (auto k = minSize; k <= maxSize; k++) {
			// Find the corresponding cluster
			comp[Species::I] = k;
			comp[Species::Loop] = 0;
			auto iProdId = subpaving.findTileId(comp);
			if (iProdId != subpaving.invalidIndex() and
				iProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, iProdId});
				if (lo1[Species::I] == 1 || lo2[Species::I] == 1) {
					this->addDissociationReaction(tag, {iProdId, i, j});
				}
				previousIndex = iProdId;

				// Find the corresponding Loop if above 40
				if (k >= 40) {
					comp[Species::I] = 0;
					comp[Species::Loop] = k;
					auto lProdId = subpaving.findTileId(comp);
					if (lProdId != subpaving.invalidIndex()) {
						// Only if both reactants are 111 loops
						if (lo1[Species::I] > 3 and lo2[Species::I] > 3)
							this->addProductionReaction(tag, {i, j, lProdId});
					}
				}
			}
		}
		return;
	}

	// V + I = V or I or 0
	if ((lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::I))) {
		// They both can be grouped
		auto vReg = lo1.isOnAxis(Species::V) ? cl1Reg : cl2Reg;
		auto iReg = lo1.isOnAxis(Species::V) ? cl2Reg : cl1Reg;

		Composition loV = vReg.getOrigin();
		Composition hiV = vReg.getUpperLimitPoint();
		Composition loI = iReg.getOrigin();
		Composition hiI = iReg.getUpperLimitPoint();

		int minSize = loV[Species::V] - hiI[Species::I] + 1;
		int maxSize = hiV[Species::V] - 1 - loI[Species::I];
		Composition comp = Composition::zero();
		// Loop on the possibilities
		for (int k = minSize; k <= maxSize; k++) {
			// 3 cases
			if (k > 0) {
				// Looking for V cluster
				comp[Species::V] = k;
				comp[Species::I] = 0;
				auto vProdId = subpaving.findTileId(comp);
				if (vProdId != subpaving.invalidIndex() and
					vProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, vProdId});
					// No dissociation
					previousIndex = vProdId;
				}
			}
			else if (k < 0) {
				// Looking for I cluster
				comp[Species::I] = -k;
				comp[Species::V] = 0;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() and
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
					previousIndex = iProdId;
				}
			}
			else {
				// No product
				this->addProductionReaction(tag, {i, j});
				// No dissociation
			}
		}
		return;
	}

	// V + V = V
	if (lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::V)) {
		// Compute the composition of the new cluster
		auto minSize = lo1[Species::V] + lo2[Species::V];
		auto maxSize = hi1[Species::V] + hi2[Species::V] - 2;

		Composition comp = Composition::zero();
		// Loop on the possibilities
		for (auto k = minSize; k <= maxSize; k++) {
			// Find the corresponding cluster
			comp[Species::V] = k;
			auto vProdId = subpaving.findTileId(comp);
			if (vProdId != subpaving.invalidIndex() and
				vProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, vProdId});
				if (lo1[Species::V] == 1 || lo2[Species::V] == 1) {
					this->addDissociationReaction(tag, {vProdId, i, j});
				}
				previousIndex = vProdId;
			}
		}
		return;
	}

	// S + S = S
	if (lo1.isOnAxis(Species::S) && lo2.isOnAxis(Species::S)) {
		// Compute the composition of the new cluster
		auto minSize = lo1[Species::S] + lo2[Species::S];
		auto maxSize = hi1[Species::S] + hi2[Species::S] - 2;

		Composition comp = Composition::zero();
		// Loop on the possibilities
		for (auto k = minSize; k <= maxSize; k++) {
			// Find the corresponding cluster
			comp[Species::S] = k;
			auto sProdId = subpaving.findTileId(comp);
			if (sProdId != subpaving.invalidIndex() and
				sProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, sProdId});
				if (lo1[Species::S] == 1 || lo2[Species::S] == 1) {
					this->addDissociationReaction(tag, {sProdId, i, j});
				}
				previousIndex = sProdId;
			}
		}
		return;
	}

	// I + Loop = Loop
	if (lo1.isOnAxis(Species::I) && lo2.isOnAxis(Species::Loop) ||
		(lo1.isOnAxis(Species::Loop) && lo2.isOnAxis(Species::I))) {
		// Compute the composition of the new cluster
		auto minSize = lo1[Species::I] + lo2[Species::I] + lo1[Species::Loop] +
			lo2[Species::Loop];
		auto maxSize = hi1[Species::I] + hi2[Species::I] + hi1[Species::Loop] +
			hi2[Species::Loop] - 4;

		Composition comp = Composition::zero();
		// Loop on the possibilities
		for (auto k = minSize; k <= maxSize; k++) {
			// Find the corresponding cluster
			comp[Species::Loop] = k;
			auto lProdId = subpaving.findTileId(comp);
			if (lProdId != subpaving.invalidIndex() and
				lProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, lProdId});
				// No Dissociation
				previousIndex = lProdId;
			}
		}
		return;
	}

	// V + Loop = Loop or I
	if ((lo1.isOnAxis(Species::Loop) && lo2.isOnAxis(Species::V)) ||
		(lo1.isOnAxis(Species::V) && lo2.isOnAxis(Species::Loop))) {
		// They both can be grouped
		auto vReg = lo1.isOnAxis(Species::V) ? cl1Reg : cl2Reg;
		auto lReg = lo1.isOnAxis(Species::V) ? cl2Reg : cl1Reg;

		Composition loV = vReg.getOrigin();
		Composition hiV = vReg.getUpperLimitPoint();
		Composition loL = lReg.getOrigin();
		Composition hiL = lReg.getUpperLimitPoint();

		int minSize = loV[Species::V] - hiL[Species::Loop] + 1;
		int maxSize = hiV[Species::V] - 1 - loL[Species::Loop];
		Composition comp = Composition::zero();
		// Loop on the possibilities
		for (int k = minSize; k <= maxSize; k++) {
			// k is always negative here
			// Looking for loop cluster
			comp[Species::Loop] = -k;
			auto lProdId = subpaving.findTileId(comp);
			if (lProdId != subpaving.invalidIndex() and
				lProdId != previousIndex) {
				this->addProductionReaction(tag, {i, j, lProdId});
				// No dissociation
				previousIndex = lProdId;
			}
			// Test if I exist
			else if (lProdId == subpaving.invalidIndex()) {
				comp[Species::Loop] = 0;
				comp[Species::I] = -k;
				auto iProdId = subpaving.findTileId(comp);
				if (iProdId != subpaving.invalidIndex() and
					iProdId != previousIndex) {
					this->addProductionReaction(tag, {i, j, iProdId});
					// No dissociation
					previousIndex = iProdId;
				}
			}
		}
		return;
	}
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
RPVReactionGenerator::addSinks(IndexType i, TTag tag) const
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
	if (clReg.isSimplex() && lo.isOnAxis(Species::V) && lo[Species::V] < 5) {
		this->addSinkReaction(tag, {i, NetworkType::invalidIndex()});
	}
}

inline ReactionCollection<RPVReactionGenerator::NetworkType>
RPVReactionGenerator::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->_clusterData.numClusters, this->_enableReadRates,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getSinkReactions());
	return ret;
}
} // namespace detail

inline detail::RPVReactionGenerator
RPVReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::RPVReactionGenerator{*this};
}
} // namespace network
} // namespace core
} // namespace xolotl
