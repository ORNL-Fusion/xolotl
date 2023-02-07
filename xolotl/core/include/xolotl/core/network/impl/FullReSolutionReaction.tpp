#pragma once

#include <xolotl/core/network/detail/ReactionUtility.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
FullReSolutionReaction<TNetwork, TDerived>::FullReSolutionReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_product(cluster1)
{
	this->asDerived()->setSize();

	this->copyMomentIds(_reactant, _reactantMomentIds);

	// static
	const auto dummyRegion = Region(Composition{});

	auto clReg = this->_clusterData->getCluster(_reactant).getRegion();
	auto prodReg = this->_clusterData->getCluster(_product).getRegion();

	_reactantVolume = clReg.volume();
	_productVolume = prodReg.volume();

	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
FullReSolutionReaction<TNetwork, TDerived>::FullReSolutionReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	FullReSolutionReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computeCoefficients()
{
	// static
	const auto dummyRegion = Region(Composition{});

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto clReg = this->_clusterData->getCluster(_reactant).getRegion();
	auto prodReg = this->_clusterData->getCluster(_product).getRegion();
	const auto& clDisp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			clReg);
	const auto& prodDisp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			prodReg);
	auto prod2Reg = dummyRegion;
	auto cl2Reg = dummyRegion;

	// Initialize the reflected regions
	auto rRegions = detail::updateReflectedRegionsForCoefs<nMomentIds>(
		prodReg, prod2Reg, clReg, cl2Reg);
	auto clRR = rRegions[2];
	auto prRR = rRegions[0];
	auto cl2RR = rRegions[3];
	auto pr2RR = rRegions[1];

	// TODO: check this is correct
	auto nOverlap = _reactantVolume;
	for (auto i : speciesRangeNoI) {
		this->_widths[i()] = clRR[i()].end() - clRR[i()].begin();
	}

	// The first coefficient is simply the overlap because it is the sum over 1
	this->_coefs(0, 0, 0, 0) = nOverlap;
	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		// First order sum
		this->_coefs(i() + 1, 0, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), clRR, cl2RR, clRR, cl2RR);
	}

	// First moments
	for (auto k : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[k()];
		// Reactant
		this->_coefs(0, 0, 0, k() + 1) =
			this->_coefs(k() + 1, 0, 0, 0) / clDisp[k()];
	}

	// Now we loop over the 1 dimension of the coefs to compute all the
	// possible sums over distances for the flux
	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		// Now we deal with the coefficients needed for the partial derivatives
		// Starting with the reactant
		for (auto k : speciesRangeNoI) {
			// Second order sum
			if (k == i) {
				this->_coefs(i() + 1, 0, 0, k() + 1) = factor *
					detail::computeSecondOrderSum(
						i(), clRR, cl2RR, clRR, cl2RR) /
					clDisp[k()];
			}
			else {
				this->_coefs(i() + 1, 0, 0, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(k() + 1, 0, 0, 0) / (nOverlap * clDisp[k()]);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
FullReSolutionReaction<TNetwork, TDerived>::computeRate(
	IndexType gridIndex, double time)
{
	const double lambda = 1.0;
	const double mu = 6000.0;

	// Compute the fraction rate
	auto cl = this->_clusterData->getCluster(_reactant);
	auto radius = cl.getReactionRadius();
	double fractionRate =
		2.0 * ::xolotl::core::pi * mu * pow(radius + lambda, 2.0);
	// Resolution rate depends on fission rate
	double resolutionRate = this->_clusterData->fissionRate();

	return fractionRate * resolutionRate;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData->getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	auto prod = this->_clusterData->getCluster(_product);
	const auto& prodReg = prod.getRegion();

	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_reactant, _reactantMomentIds[i()], connectivity);
			this->addConnectivity(
				_reactantMomentIds[i()], _reactant, connectivity);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					this->addConnectivity(_reactantMomentIds[i()],
						_reactantMomentIds[j()], connectivity);
				}
			}
		}
	}
	// Product connects with  the reactant
	this->addConnectivity(_product, _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_product, _reactantMomentIds[i()], connectivity);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData->getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	auto prod = this->_clusterData->getCluster(_product);
	const auto& prodReg = prod.getRegion();

	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (i() == j())
					this->addConnectivity(_reactantMomentIds[i()],
						_reactantMomentIds[j()], connectivity);
			}
		}
	}
	// Product connects with  the reactant
	if (_product == _reactant)
		this->addConnectivity(_product, _reactant, connectivity);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Initialize the concentrations that will be used in the loops
	auto cR = concentrations[_reactant];
	Kokkos::Array<double, nMomentIds> cmR;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] == invalidIndex) {
			cmR[i()] = 0.0;
		}
		else
			cmR[i()] = concentrations[_reactantMomentIds[i()]];
	}

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * cR;
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) * cmR[i()];
	}
	f *= this->_rate(gridIndex);
	Kokkos::atomic_sub(&fluxes[_reactant], f / _reactantVolume);
	Kokkos::atomic_add(&fluxes[_product], _size * f / _productVolume);

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		// First for the reactant
		if (_reactantMomentIds[k()] != invalidIndex) {
			f = this->_coefs(0, 0, 0, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 0, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[k()]], f / _reactantVolume);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	using AmountType = typename NetworkType::AmountType;
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = this->_rate(gridIndex) / _reactantVolume;
	// Compute the values
	Kokkos::atomic_sub(
		&values(_connEntries[0][0][0][0]), df * this->_coefs(0, 0, 0, 0));
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_sub(&values(_connEntries[0][0][0][1 + i()]),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}

	// For the product
	df = this->_rate(gridIndex) / _productVolume;
	Kokkos::atomic_add(&values(_connEntries[1][0][0][0]),
		_size * df * this->_coefs(0, 0, 0, 0));
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_add(&values(_connEntries[1][0][0][1 + i()]),
				_size * df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = this->_rate(gridIndex) / _reactantVolume;
			// Compute the values
			Kokkos::atomic_sub(&values(_connEntries[0][1 + k()][0][0]),
				df * this->_coefs(0, 0, 0, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 0, k() + 1));
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	using AmountType = typename NetworkType::AmountType;
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = this->_rate(gridIndex) / _reactantVolume;
	// Compute the values
	Kokkos::atomic_sub(
		&values(_connEntries[0][0][0][0]), df * this->_coefs(0, 0, 0, 0));
	// For the first product
	df = this->_rate(gridIndex) / _productVolume;
	if (_product == _reactant)
		Kokkos::atomic_add(&values(_connEntries[1][0][0][0]),
			_size * df * this->_coefs(0, 0, 0, 0));

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = this->_rate(gridIndex) / _reactantVolume;
			// Compute the values
			for (auto i : speciesRangeNoI) {
				if (k() == i())
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 0, k() + 1));
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::computeConstantRates(
	ConcentrationsView concentrations, RatesView rates, BelongingView isInSub,
	OwnedSubMapView backMap, IndexType gridIndex)
{
	return;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
FullReSolutionReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// This type of reaction doesn't count
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
FullReSolutionReaction<TNetwork, TDerived>::mapJacobianEntries(
	Connectivity connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	_connEntries[0][0][0][0] = connectivity(_reactant, _reactant);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			_connEntries[0][0][0][1 + i()] =
				connectivity(_reactant, _reactantMomentIds[i()]);
			_connEntries[0][1 + i()][0][0] =
				connectivity(_reactantMomentIds[i()], _reactant);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					_connEntries[0][1 + i()][0][1 + j()] = connectivity(
						_reactantMomentIds[i()], _reactantMomentIds[j()]);
				}
			}
		}
	}

	_connEntries[1][0][0][0] = connectivity(_product, _reactant);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			_connEntries[1][0][0][1 + i()] =
				connectivity(_product, _reactantMomentIds[i()]);
		}
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
