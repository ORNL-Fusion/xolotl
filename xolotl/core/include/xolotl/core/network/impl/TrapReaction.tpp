#pragma once

#include <plsm/EnumIndexed.h>

#include <xolotl/core/Constants.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION void
TrapReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	// Get the concentrations for mobile and trapped
	auto cm = concentrations[_reactant];
	auto ct = concentrations[_trapped];

	// Get the rates and other parameters
	auto trappingRate = this->computeTrappingRate(gridIndex);
	auto deTrappingRate = this->computeDeTrappingRate(gridIndex);
	auto trapDensity = this->asDerived()->getDensity();

	// Compute the total trapping flux
	double flux = trappingRate * (trapDensity - ct) * cm;
	flux -= deTrappingRate * ct;

	// The trapped concentration increases
	Kokkos::atomic_add(&fluxes[_trapped], flux);
	// The mobile concentration decreases
	Kokkos::atomic_sub(&fluxes[_reactant], flux);
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION void
TrapReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Get the concentrations for mobile and trapped
	auto cm = concentrations[_reactant];
	auto ct = concentrations[_trapped];

	// Get the rates and other parameters
	auto trappingRate = this->computeTrappingRate(gridIndex);
	auto deTrappingRate = this->computeDeTrappingRate(gridIndex);
	auto trapDensity = this->asDerived()->getDensity();

	// Compute the partial with respect to the mobile concentration
	double partial = trappingRate * (trapDensity - ct);
	Kokkos::atomic_sub(&values(_connEntries[0][0]), partial);
	Kokkos::atomic_add(&values(_connEntries[1][0]), partial);

	// Compute the partial with respect to the trapped concentration
	partial = trappingRate * cm + deTrappingRate;
	Kokkos::atomic_add(&values(_connEntries[0][1]), partial);
	Kokkos::atomic_sub(&values(_connEntries[1][1]), partial);
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION void
TrapReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Get the concentrations for mobile and trapped
	auto cm = concentrations[_reactant];
	auto ct = concentrations[_trapped];

	// Get the rates and other parameters
	auto trappingRate = this->computeTrappingRate(gridIndex);
	auto deTrappingRate = this->computeDeTrappingRate(gridIndex);
	auto trapDensity = this->asDerived()->getDensity();

	// Compute the partial with respect to the mobile concentration
	double partial = trappingRate * (trapDensity - ct);
	Kokkos::atomic_sub(&values(_connEntries[0][0]), partial);

	// Compute the partial with respect to the trapped concentration
	partial = trappingRate * cm + deTrappingRate;
	Kokkos::atomic_sub(&values(_connEntries[1][1]), partial);
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION double
TrapReaction<TNetwork, TDerived>::computeTrappingRate(
	IndexType gridIndex, double time)
{
	auto cl = this->_clusterData->getCluster(_reactant);
	double dc = cl.getDiffusionCoefficient(gridIndex);

	double rate = this->asDerived()->getStrength() * dc;

	return rate;
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION double
TrapReaction<TNetwork, TDerived>::computeDeTrappingRate(
	IndexType gridIndex, double time)
{
	double temperature = this->_clusterData->temperature(gridIndex);
	constexpr double k_B = ::xolotl::core::kBoltzmann;

	double rate = this->asDerived()->getFrequency();
	rate *= exp(-this->asDerived()->getEnergy() / (temperature * k_B));

	return rate;
}
} // namespace network
} // namespace core
} // namespace xolotl
