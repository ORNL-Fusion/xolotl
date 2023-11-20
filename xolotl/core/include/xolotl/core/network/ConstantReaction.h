#pragma once

#include <xolotl/core/network/Reaction.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>

namespace xolotl
{
namespace core
{
namespace network
{
/**
 * @brief Class implementing reaction where the other cluster is missing from
 * the network but its rate is accounted for as a constant rate that can be
 * updated.
 *
 * @tparam TNetwork The network type
 * @tparam TDerived The derived class type.
 */
template <typename TNetwork, typename TDerived>
class ConstantReaction : public Reaction<TNetwork, TDerived>
{
	friend class Reaction<TNetwork, TDerived>;

public:
	using NetworkType = TNetwork;
	using Superclass = Reaction<TNetwork, TDerived>;
	using IndexType = typename Superclass::IndexType;
	using Connectivity = typename Superclass::Connectivity;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;
	using RatesView = typename Superclass::RatesView;
	using ConnectivitiesView = typename Superclass::ConnectivitiesView;
	using ConnectivitiesPairView = typename Superclass::ConnectivitiesPairView;
	using BelongingView = typename Superclass::BelongingView;
	using OwnedSubMapView = typename Superclass::OwnedSubMapView;
	using AmountType = typename Superclass::AmountType;
	using ReactionDataRef = typename Superclass::ReactionDataRef;
	using ClusterData = typename Superclass::ClusterData;
	using RateVector = IReactionNetwork::RateVector;

	ConstantReaction() = default;

	KOKKOS_INLINE_FUNCTION
	ConstantReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		IndexType cluster0, IndexType cluster1) :
		Superclass(reactionData, clusterData, reactionId),
		_reactants({cluster0, cluster1})
	{
		for (auto i : {0, 1}) {
			this->copyMomentIds(_reactants[i], _reactantMomentIds[i]);
		}
		this->initialize();
	}

	KOKKOS_INLINE_FUNCTION
	ConstantReaction(ReactionDataRef reactionData,
		const ClusterData& clusterData, IndexType reactionId,
		const detail::ClusterSet& clusterSet) :
		ConstantReaction(reactionData, clusterData, reactionId,
			clusterSet.cluster0, clusterSet.cluster1)
	{
	}

	static detail::CoefficientsView
	allocateCoefficientsView(IndexType)
	{
		return detail::CoefficientsView();
	}

	static detail::ConstantRateView
	allocateConstantRateView(IndexType size, IndexType gridSize)
	{
		return detail::ConstantRateView("Constant Rates", size, gridSize,
			Superclass::coeffsSingleExtent, Superclass::coeffsSingleExtent);
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeRate(IndexType gridIndex, double time = 0.0)
	{
		return 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void
	setRate(RatesView rates, IndexType gridIndex)
	{
		auto dof = rates.extent(0);
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		if (_reactants[1] == invalidIndex) {
			this->_constantRates(gridIndex, 0, 0) = rates(_constantRateEntries[0][0]);
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					this->_constantRates(gridIndex, 1 + i(), 0) =
						rates(_constantRateEntries[1 + i()][0]);
				}
			}
		}
		else {
			this->_constantRates(gridIndex, 0, 0) = rates(_constantRateEntries[0][0]);
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[1][i()] != invalidIndex) {
					this->_constantRates(gridIndex, 0, 1 + i()) =
						rates(_constantRateEntries[0][1 + i()]);
				}
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					this->_constantRates(gridIndex, 1 + i(), 0) =
						rates(_constantRateEntries[1 + i()][0]);
					for (auto j : speciesRangeNoI) {
						if (_reactantMomentIds[1][j()] != invalidIndex) {
							this->_constantRates(gridIndex, 1 + i(), 1 + j()) =
								rates(_constantRateEntries[1 + i()][1 + j()]);
						}
					}
				}
			}
		}
	}

private:
	KOKKOS_INLINE_FUNCTION
	void
	computeCoefficients()
	{
		// No coefs
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeConnectivity(const Connectivity& connectivity)
	{
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		// Reactant 1 with reactant 1
		this->addConnectivity(_reactants[0], _reactants[0], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				this->addConnectivity(
					_reactants[0], _reactantMomentIds[0][i()], connectivity);
				this->addConnectivity(
					_reactantMomentIds[0][i()], _reactants[0], connectivity);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[0][j()] != invalidIndex) {
						this->addConnectivity(_reactantMomentIds[0][i()],
							_reactantMomentIds[0][j()], connectivity);
					}
				}
			}
		}

		if (_reactants[1] == invalidIndex)
			return;

		// Reactant 2 with reactant 1
		this->addConnectivity(_reactants[0], _reactants[1], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] != invalidIndex) {
				this->addConnectivity(
					_reactants[0], _reactantMomentIds[1][i()], connectivity);
			}
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				this->addConnectivity(
					_reactantMomentIds[0][i()], _reactants[1], connectivity);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[1][j()] != invalidIndex) {
						this->addConnectivity(_reactantMomentIds[0][i()],
							_reactantMomentIds[1][j()], connectivity);
					}
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedConnectivity(const Connectivity& connectivity)
	{
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		// Each reactant connects with all the reactants
		// Reactant 1 with reactant 1
		this->addConnectivity(_reactants[0], _reactants[0], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				for (auto j : speciesRangeNoI) {
					if (i() == j())
						this->addConnectivity(_reactantMomentIds[0][i()],
							_reactantMomentIds[0][j()], connectivity);
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeFlux(ConcentrationsView concentrations, FluxesView fluxes,
		IndexType gridIndex)
	{
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		if (_reactants[1] == invalidIndex) {
			Kokkos::atomic_add(
				&fluxes(_reactants[0]), this->_constantRates(gridIndex, 0, 0));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					Kokkos::atomic_add(&fluxes(_reactantMomentIds[0][i()]),
						this->_constantRates(gridIndex, 1 + i(), 0));
				}
			}
		}
		else {
			Kokkos::atomic_add(&fluxes(_reactants[0]),
				concentrations(_reactants[1]) *
					this->_constantRates(gridIndex, 0, 0));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[1][i()] != invalidIndex) {
					Kokkos::atomic_add(&fluxes(_reactants[0]),
						concentrations(_reactantMomentIds[1][i()]) *
							this->_constantRates(gridIndex, 0, 1 + i()));
				}
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					Kokkos::atomic_add(&fluxes(_reactantMomentIds[0][i()]),
						concentrations(_reactants[1]) *
							this->_constantRates(gridIndex, 1 + i(), 0));
					for (auto j : speciesRangeNoI) {
						if (_reactantMomentIds[1][j()] != invalidIndex) {
							Kokkos::atomic_add(
								&fluxes(_reactantMomentIds[0][i()]),
								concentrations(_reactantMomentIds[1][j()]) *
									this->_constantRates(
										gridIndex, 1 + i(), 1 + j()));
						}
					}
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	computePartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		if (_reactants[1] == invalidIndex)
			return;

		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		Kokkos::atomic_add(&values(_connEntries[0][0][0][0]),
			this->_constantRates(gridIndex, 0, 0));
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] != invalidIndex) {
				Kokkos::atomic_add(&values(_connEntries[0][0][0][1 + i()]),
					this->_constantRates(gridIndex, 0, 1 + i()));
			}
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				Kokkos::atomic_add(&values(_connEntries[0][1 + i()][0][0]),
					this->_constantRates(gridIndex, 1 + i(), 0));
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[1][j()] != invalidIndex) {
						Kokkos::atomic_add(
							&values(_connEntries[0][1 + i()][0][1 + j()]),
							this->_constantRates(gridIndex, 1 + i(), 1 + j()));
					}
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeReducedPartialDerivatives(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex)
	{
		if (_reactants[1] == invalidIndex)
			return;

		if (_reactants[1] != _reactants[0])
			return;

		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		Kokkos::atomic_add(&values(_connEntries[0][0][0][0]),
			this->_constantRates(gridIndex, 0, 0));
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[1][j()] != invalidIndex) {
					Kokkos::atomic_add(
						&values(_connEntries[0][1 + i()][0][1 + j()]),
						this->_constantRates(gridIndex, 1 + i(), 1 + j()));
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	computeConstantRates(ConcentrationsView concentrations, RatesView rates,
		BelongingView isInSub, IndexType subId, IndexType gridIndex)
	{
		return;
	}

	KOKKOS_INLINE_FUNCTION
	void
	getConstantConnectivities(ConnectivitiesView conns, BelongingView isInSub,
		OwnedSubMapView backMap)
	{
		return;
	}

	KOKKOS_INLINE_FUNCTION
	double
	computeLeftSideRate(ConcentrationsView concentrations, IndexType clusterId,
		IndexType gridIndex)
	{
		return 0.0;
	}

	KOKKOS_INLINE_FUNCTION
	void
	mapJacobianEntries(Connectivity connectivity)
	{
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		if (_reactants[1] == invalidIndex)
			return;

		_connEntries[0][0][0][0] = connectivity(_reactants[0], _reactants[1]);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] != invalidIndex) {
				_connEntries[0][0][0][1 + i()] =
					connectivity(_reactants[0], _reactantMomentIds[1][i()]);
			}
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				_connEntries[0][1 + i()][0][0] =
					connectivity(_reactantMomentIds[0][i()], _reactants[1]);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[1][j()] != invalidIndex) {
						_connEntries[0][1 + i()][0][1 + j()] =
							connectivity(_reactantMomentIds[0][i()],
								_reactantMomentIds[1][j()]);
					}
				}
			}
		}
	}

	KOKKOS_INLINE_FUNCTION
	void
	mapRateEntries(ConnectivitiesPairView connectivityRow,
		ConnectivitiesPairView connectivityEntries, BelongingView isInSub,
		OwnedSubMapView backMap, IndexType subId)
	{
		printf("constant");
		auto dof = connectivityRow.extent(0) - 1;
		constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

		if (_reactants[1] == invalidIndex) {
			_constantRateEntries[0][0] = this->getPosition(
				_reactants[0], dof, connectivityRow, connectivityEntries);
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					_constantRateEntries[1 + i()][0] =
						this->getPosition(_reactantMomentIds[0][i()], dof,
							connectivityRow, connectivityEntries);
				}
			}
		}
		else {
			_constantRateEntries[0][0] = this->getPosition(_reactants[0], _reactants[1],
				connectivityRow, connectivityEntries);
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[1][i()] != invalidIndex) {
					_constantRateEntries[0][1 + i()] = this->getPosition(_reactants[0],
						_reactantMomentIds[1][i()], connectivityRow,
						connectivityEntries);
				}
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					_constantRateEntries[1 + i()][0] = this->getPosition(
						_reactantMomentIds[0][i()], _reactants[1],
						connectivityRow, connectivityEntries);
					for (auto j : speciesRangeNoI) {
						if (_reactantMomentIds[1][j()] != invalidIndex) {
							_constantRateEntries[1 + i()][1 + j()] =
								this->getPosition(_reactantMomentIds[0][i()],
									_reactantMomentIds[1][j()], connectivityRow,
									connectivityEntries);
						}
					}
				}
			}
		}
	}

protected:
	static constexpr auto invalidIndex = Superclass::invalidIndex;
	util::Array<IndexType, 2> _reactants{invalidIndex, invalidIndex};

	static constexpr auto nMomentIds = Superclass::nMomentIds;
	util::Array<IndexType, 2, nMomentIds> _reactantMomentIds;

	util::Array<IndexType, 1, 1 + nMomentIds, 1, 1 + nMomentIds> _connEntries;
	util::Array<IndexType, Superclass::coeffsSingleExtent,
		Superclass::coeffsSingleExtent>
		_constantRateEntries;
};
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/detail/ConstantReactionGenerator.h>
