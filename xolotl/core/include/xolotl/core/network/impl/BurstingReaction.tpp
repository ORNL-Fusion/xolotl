#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
BurstingReaction<TNetwork, TDerived>::BurstingReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_product(cluster1)
{
	this->copyMomentIds(_reactant, _reactantMomentIds);

	auto clReg = this->_clusterData->getCluster(_reactant).getRegion();
	_reactantVolume = clReg.volume();

	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
BurstingReaction<TNetwork, TDerived>::BurstingReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	BurstingReaction(reactionData, clusterData, reactionId, clusterSet.cluster0,
		clusterSet.cluster1)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeCoefficients()
{
	using Species = typename TNetwork::Species;
	// static
	const auto dummyRegion = Region(Composition{});

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto clReg = this->_clusterData->getCluster(_reactant).getRegion();
	const auto& prReg = (_product == invalidIndex) ?
		dummyRegion :
		this->_clusterData->getCluster(_product).getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();
	Composition prod = prReg.getOrigin();
	auto nV = prod[Species::V];

	this->_coefs(0, 0, 0, 0) =
		clReg.volume() / (hi[Species::V] - lo[Species::V]);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			if (isVacancy(i)) {
				auto factor = this->_coefs(0, 0, 0, 0) *
					(nV - (hi[Species::V] + lo[Species::V] - 1) / 2.0);
				this->_coefs(i() + 1, 0, 0, 0) = factor;

				auto disp = clReg.dispersion();
				factor = this->_coefs(0, 0, 0, 0) *
					(nV - (hi[Species::V] + lo[Species::V] - 1) / 2.0) /
					disp[i()];
				this->_coefs(0, 0, 0, i() + 1) = factor;

				factor = this->_coefs(0, 0, 0, 0) *
					(nV - (hi[Species::V] + lo[Species::V] - 1) / 2.0) *
					(nV - (hi[Species::V] + lo[Species::V] - 1) / 2.0) /
					disp[i()];
				this->_coefs(i() + 1, 0, 0, i() + 1) = factor;
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
BurstingReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::updateRates(double largestRate)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

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
	if (_product != invalidIndex) {
		this->addConnectivity(_product, _reactant, connectivity);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] != invalidIndex) {
				this->addConnectivity(
					_product, _reactantMomentIds[i()], connectivity);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

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
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeFlux(
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

	auto rate = this->asDerived()->getAppliedRate(gridIndex);

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * cR;
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) * cmR[i()];
	}
	f *= rate;

	Kokkos::atomic_sub(&fluxes[_reactant], f / (double)_reactantVolume);
	if (_product != invalidIndex)
		Kokkos::atomic_add(&fluxes[_product], f);

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		// For the reactant
		if (_reactantMomentIds[k()] != invalidIndex) {
			f = this->_coefs(0, 0, 0, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 0, k() + 1) * cmR[i()];
			}
			f *= rate;
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[k()]], f / (double)_reactantVolume);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto rate = this->asDerived()->getAppliedRate(gridIndex);

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = rate / (double)_reactantVolume;
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
	if (_product != invalidIndex) {
		df = rate;
		Kokkos::atomic_add(
			&values(_connEntries[1][0][0][0]), df * this->_coefs(0, 0, 0, 0));

		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] != invalidIndex) {
				Kokkos::atomic_add(&values(_connEntries[1][0][0][1 + i()]),
					df * this->_coefs(i() + 1, 0, 0, 0));
			}
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = rate / (double)_reactantVolume;
			// Compute the values
			Kokkos::atomic_sub(&values(_connEntries[0][1 + k()][0][0]),
				df * this->_coefs(0, 0, 0, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						df * this->_coefs(k() + 1, 0, 0, i() + 1));
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto rate = this->asDerived()->getAppliedRate(gridIndex);
	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]),
		rate * this->_coefs(0, 0, 0, 0) / (double)_reactantVolume);

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			auto df = rate / (double)_reactantVolume;
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] == _reactantMomentIds[k()]) {
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						df * this->_coefs(k() + 1, 0, 0, i() + 1));
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
BurstingReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::mapJacobianEntries(
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

	if (_product != invalidIndex) {
		_connEntries[1][0][0][0] = connectivity(_product, _reactant);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] != invalidIndex) {
				_connEntries[1][0][0][1 + i()] =
					connectivity(_product, _reactantMomentIds[i()]);
			}
		}
	}
}

} // namespace network
} // namespace core
} // namespace xolotl
