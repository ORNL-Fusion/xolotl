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
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_product(cluster1)
{
	this->copyMomentIds(_reactant, _reactantMomentIds);
	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
BurstingReaction<TNetwork, TDerived>::BurstingReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
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

	auto clReg = this->_clusterData.getCluster(_reactant).getRegion();
	const auto& prReg = (_product == invalidIndex) ?
		dummyRegion :
		this->_clusterData.getCluster(_product).getRegion();
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
	return 1.0e9;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::updateRates(double largestRate)
{
	for (IndexType i = 0; i < this->_rate.extent(0); ++i) {
		this->_rate(i) = 1.0e9;
	}
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
double
BurstingReaction<TNetwork, TDerived>::getAppliedRate(IndexType gridIndex) const
{
	// Get the radius of the cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	auto radius = cl.getReactionRadius();

	// Get the current depth
	auto depth = this->_clusterData.getDepth();
	auto tau = this->_clusterData.getTauBursting();
	auto f = this->_clusterData.getFBursting();
	return f * (radius / depth) *
		util::min(1.0, exp(-(depth - tau) / (2.0 * tau)));
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	AmountType volCl = clReg.volume();
	AmountType volProd = 1;

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

	auto rate = getAppliedRate(gridIndex);

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * cR;
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) * cmR[i()];
	}
	f *= rate;

	Kokkos::atomic_sub(&fluxes[_reactant], f / (double)volCl);
	if (_product != invalidIndex)
		Kokkos::atomic_add(&fluxes[_product], f / (double)volProd);

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
				&fluxes[_reactantMomentIds[k()]], f / (double)volCl);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
BurstingReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	AmountType volCl = clReg.volume();
	AmountType volProd = 1;

	auto rate = getAppliedRate(gridIndex);

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = rate / (double)volCl;
	// Compute the values
	Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)),
		df * this->_coefs(0, 0, 0, 0));
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_sub(
				&values(connectivity(_reactant, _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the product
	if (_product != invalidIndex) {
		df = rate / (double)volProd;
		Kokkos::atomic_add(&values(connectivity(_product, _reactant)),
			df * this->_coefs(0, 0, 0, 0));

		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] != invalidIndex) {
				Kokkos::atomic_add(
					&values(connectivity(_product, _reactantMomentIds[i()])),
					df * this->_coefs(i() + 1, 0, 0, 0));
			}
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = rate / (double)volCl;
			// Compute the values
			Kokkos::atomic_sub(
				&values(connectivity(_reactantMomentIds[k()], _reactant)),
				df * this->_coefs(0, 0, 0, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_sub(
						&values(connectivity(
							_reactantMomentIds[k()], _reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 0, k() + 1));
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
	Connectivity connectivity, IndexType gridIndex)
{
	auto rate = getAppliedRate(gridIndex);
	Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)), rate);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
BurstingReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
