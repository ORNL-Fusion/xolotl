#pragma once

#include <xolotl/core/network/impl/NucleationReaction.tpp>
#include <xolotl/core/network/impl/ReSolutionReaction.tpp>
#include <xolotl/core/network/impl/Reaction.tpp>
#include <xolotl/core/network/impl/SinkReaction.tpp>

namespace xolotl
{
namespace core
{
namespace network
{
namespace ne
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
double
getRate(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg, const double r0,
	const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;

	double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1) * ::xolotl::core::zFactor;

	return kPlus;
}
} // namespace ne

KOKKOS_INLINE_FUNCTION
double
NEProductionReaction::getRateForProduction(IndexType gridIndex)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionFactor();
	double dc1 = cl1.getDiffusionFactor();

	return ne::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
NEProductionReaction::computeRate(IndexType gridIndex, double time)
{
	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	Composition cl0Comp = cl0.getRegion().getOrigin();
	Composition cl1Comp = cl1.getRegion().getOrigin();

	//		if (cl0Comp[Species::Xe] > 0 and cl1Comp[Species::Xe] > 0)
	//			return 0.0;

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	// Read the rates if available
	double rate = -1.0;
	if (this->_clusterData->extraData.fileClusterMap.exists(_reactants[0]) and
		this->_clusterData->extraData.fileClusterMap.exists(_reactants[1])) {
		auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
		rate = this->_clusterData->extraData.constantRates(
			clusterMap.value_at(clusterMap.find(_reactants[0])),
			clusterMap.value_at(clusterMap.find(_reactants[1])), 0);

		if (this->_deltaG0 > 0.0) {
			rate = this->_clusterData->extraData.constantRates(
					   clusterMap.value_at(clusterMap.find(_reactants[0])),
					   clusterMap.value_at(clusterMap.find(_reactants[1])), 1) *
				exp(-this->_deltaG0 /
					(::xolotl::core::kBoltzmann *
						this->_clusterData->temperature(gridIndex)));
		}
	}
	if (rate > 0.0) {
		return rate * 4.0 * ::xolotl::core::pi * r0 * ::xolotl::core::zFactor;
	}

	rate = getRateForProduction(gridIndex);

	return rate;
}

KOKKOS_INLINE_FUNCTION
void
NEProductionReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
	}

	if (nProd == 0) {
		double omega = this->_clusterData->atomicVolume();

		auto cl0 = this->_clusterData->getCluster(_reactants[0]);
		auto cl1 = this->_clusterData->getCluster(_reactants[1]);

		double r0 = cl0.getReactionRadius();
		double r1 = cl1.getReactionRadius();

		// Compute the flux for the 0th order moments
		auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
		double f = (1.0 / (omega * omega)) *
			this->_clusterData->extraData.constantRates(
				clusterMap.value_at(clusterMap.find(_reactants[0])),
				clusterMap.value_at(clusterMap.find(_reactants[1])), 0) *
			exp(this->_deltaG0 / ::xolotl::core::kBoltzmann *
				this->_clusterData->temperature(gridIndex)) *
			4.0 * ::xolotl::core::pi * r0 * ::xolotl::core::zFactor;

		Kokkos::atomic_add(
			&fluxes[_reactants[0]], f / (double)_reactantVolumes[0]);
		Kokkos::atomic_add(
			&fluxes[_reactants[1]], f / (double)_reactantVolumes[1]);
	}

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Initialize the concentrations that will be used in the loops
	auto cR1 = concentrations[_reactants[0]];
	Kokkos::Array<double, nMomentIds> cmR1;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[0][i()] == invalidIndex) {
			cmR1[i()] = 0.0;
		}
		else
			cmR1[i()] = concentrations[_reactantMomentIds[0][i()]];
	}
	auto cR2 = concentrations[_reactants[1]];
	Kokkos::Array<double, nMomentIds> cmR2;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] == invalidIndex) {
			cmR2[i()] = 0.0;
		}
		else
			cmR2[i()] = concentrations[_reactantMomentIds[1][i()]];
	}

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * cR1 * cR2;
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()] * cR2;
		f += this->_coefs(0, i() + 1, 0, 0) * cR1 * cmR2[i()];
		for (auto j : speciesRangeNoI) {
			f += this->_coefs(i() + 1, j() + 1, 0, 0) * cmR1[i()] * cmR2[j()];
		}
	}
	f *= this->_rate(gridIndex);

	Kokkos::atomic_sub(&fluxes[_reactants[0]], f / (double)_reactantVolumes[0]);
	Kokkos::atomic_sub(&fluxes[_reactants[1]], f / (double)_reactantVolumes[1]);
	IndexType p = 0;
	for (auto prodId : _products) {
		if (prodId == invalidIndex) {
			continue;
		}
		Kokkos::atomic_add(&fluxes[prodId], f / (double)_productVolumes[p]);
		p++;
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		// First for the first reactant
		if (_reactantMomentIds[0][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 0, k() + 1) * cR1 * cR2;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 0, k() + 1) * cmR1[i()] * cR2;
				f += this->_coefs(0, i() + 1, 0, k() + 1) * cR1 * cmR2[i()];
				for (auto j : speciesRangeNoI) {
					f += this->_coefs(i() + 1, j() + 1, 0, k() + 1) *
						cmR1[i()] * cmR2[j()];
				}
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(&fluxes[_reactantMomentIds[0][k()]],
				f / (double)_reactantVolumes[0]);
		}

		// For the second reactant
		if (_reactantMomentIds[1][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 1, k() + 1) * cR1 * cR2;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 1, k() + 1) * cmR1[i()] * cR2;
				f += this->_coefs(0, i() + 1, 1, k() + 1) * cR1 * cmR2[i()];
				for (auto j : speciesRangeNoI) {
					f += this->_coefs(i() + 1, j() + 1, 1, k() + 1) *
						cmR1[i()] * cmR2[j()];
				}
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(&fluxes[_reactantMomentIds[1][k()]],
				f / (double)_reactantVolumes[1]);
		}

		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}

			if (_productMomentIds[p][k()] != invalidIndex) {
				f = this->_coefs(0, 0, p + 2, k() + 1) * cR1 * cR2;
				for (auto i : speciesRangeNoI) {
					f += this->_coefs(i() + 1, 0, p + 2, k() + 1) * cmR1[i()] *
						cR2;
					f += this->_coefs(0, i() + 1, p + 2, k() + 1) * cR1 *
						cmR2[i()];
					for (auto j : speciesRangeNoI) {
						f += this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) *
							cmR1[i()] * cmR2[j()];
					}
				}
				f *= this->_rate(gridIndex);
				Kokkos::atomic_add(&fluxes[_productMomentIds[p][k()]],
					f / (double)_productVolumes[p]);
			}
		}
	}
}

KOKKOS_INLINE_FUNCTION
void
NEProductionReaction::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Initialize the concentrations that will be used in the loops
	auto cR1 = concentrations[_reactants[0]];
	Kokkos::Array<double, nMomentIds> cmR1;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[0][i()] == invalidIndex) {
			cmR1[i()] = 0.0;
		}
		else
			cmR1[i()] = concentrations[_reactantMomentIds[0][i()]];
	}
	auto cR2 = concentrations[_reactants[1]];
	Kokkos::Array<double, nMomentIds> cmR2;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] == invalidIndex) {
			cmR2[i()] = 0.0;
		}
		else
			cmR2[i()] = concentrations[_reactantMomentIds[1][i()]];
	}

	// Compute the partials for the 0th order moments
	// Compute the values (d / dL_0^A)
	double temp = this->_coefs(0, 0, 0, 0) * cR2;
	for (auto i : speciesRangeNoI) {
		temp += this->_coefs(0, i() + 1, 0, 0) * cmR2[i()];
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]),
		this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);
	// Second reactant
	Kokkos::atomic_sub(&values(_connEntries[1][0][0][0]),
		this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][0][0]),
			this->_rate(gridIndex) * temp / (double)_productVolumes[p]);
	}

	// Compute the values (d / dL_0^B)
	temp = this->_coefs(0, 0, 0, 0) * cR1;
	for (auto i : speciesRangeNoI) {
		temp += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()];
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(_connEntries[0][0][1][0]),
		this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);
	// Second reactant
	Kokkos::atomic_sub(&values(_connEntries[1][0][1][0]),
		this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][1][0]),
			this->_rate(gridIndex) * temp / (double)_productVolumes[p]);
	}

	for (auto i : speciesRangeNoI) {
		// (d / dL_1^A)
		if (_reactantMomentIds[0][i()] != invalidIndex) {
			temp = this->_coefs(i() + 1, 0, 0, 0) * cR2;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(i() + 1, j() + 1, 0, 0) * cmR2[j()];
			}
			// First reactant
			Kokkos::atomic_sub(&values(_connEntries[0][0][0][1 + i()]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);
			// second reactant
			Kokkos::atomic_sub(&values(_connEntries[1][0][0][1 + i()]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][0][0][1 + i()]),
					this->_rate(gridIndex) * temp / (double)_productVolumes[p]);
			}
		}

		// (d / dL_1^B)
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			temp = this->_coefs(0, i() + 1, 0, 0) * cR1;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, i() + 1, 0, 0) * cmR1[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[0][0][1][1 + i()]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);
			Kokkos::atomic_sub(&values(_connEntries[1][0][1][1 + i()]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][0][1][1 + i()]),
					this->_rate(gridIndex) * temp / (double)_productVolumes[p]);
			}
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[0][k()] != invalidIndex) {
			// First for the first reactant
			// (d / dL_0^A)
			temp = this->_coefs(0, 0, 0, k() + 1) * cR2;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(0, j() + 1, 0, k() + 1) * cmR2[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[0][1 + k()][0][0]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);

			// (d / dL_0^B)
			temp = this->_coefs(0, 0, 0, k() + 1) * cR1;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, 0, 0, k() + 1) * cmR1[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[0][1 + k()][1][0]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);

			for (auto i : speciesRangeNoI) {
				// (d / dL_1^A)
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					temp = this->_coefs(i() + 1, 0, 0, k() + 1) * cR2;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, 0, k() + 1) *
							cmR2[j()];
					}
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[0]);
				}

				// (d / dL_1^B)
				if (_reactantMomentIds[1][i()] != invalidIndex) {
					temp = this->_coefs(0, i() + 1, 0, k() + 1) * cR1;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, 0, k() + 1) *
							cmR1[j()];
					}
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][1][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[0]);
				}
			}
		}

		if (_reactantMomentIds[1][k()] != invalidIndex) {
			// First for the second reactant
			// (d / dL_0^A)
			temp = this->_coefs(0, 0, 1, k() + 1) * cR2;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(0, j() + 1, 1, k() + 1) * cmR2[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[1][1 + k()][0][0]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);

			// (d / dL_0^B)
			temp = this->_coefs(0, 0, 1, k() + 1) * cR1;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, 0, 1, k() + 1) * cmR1[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[1][1 + k()][1][0]),
				this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);

			for (auto i : speciesRangeNoI) {
				// (d / dL_1^A)
				if (_reactantMomentIds[0][i()] != invalidIndex) {
					temp = this->_coefs(i() + 1, 0, 1, k() + 1) * cR2;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, 1, k() + 1) *
							cmR2[j()];
					}
					Kokkos::atomic_sub(
						&values(_connEntries[1][1 + k()][0][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[1]);
				}

				// (d / dL_1^B)
				if (_reactantMomentIds[1][i()] != invalidIndex) {
					temp = this->_coefs(0, i() + 1, 1, k() + 1) * cR1;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, 1, k() + 1) *
							cmR1[j()];
					}
					Kokkos::atomic_sub(
						&values(_connEntries[1][1 + k()][1][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[1]);
				}
			}
		}
	}

	// Loop on the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			if (_productMomentIds[p][k()] != invalidIndex) {
				// (d / dL_0^A)
				temp = this->_coefs(0, 0, p + 2, k() + 1) * cR2;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(0, j() + 1, p + 2, k() + 1) * cmR2[j()];
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][1 + k()][0][0]),
					this->_rate(gridIndex) * temp / (double)_productVolumes[p]);

				// (d / dL_0^B)
				temp = this->_coefs(0, 0, p + 2, k() + 1) * cR1;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(j() + 1, 0, p + 2, k() + 1) * cmR1[j()];
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][1 + k()][1][0]),
					this->_rate(gridIndex) * temp / (double)_productVolumes[p]);

				for (auto i : speciesRangeNoI) {
					// (d / dL_1^A)
					if (_reactantMomentIds[0][i()] != invalidIndex) {
						temp = this->_coefs(i() + 1, 0, p + 2, k() + 1) * cR2;
						for (auto j : speciesRangeNoI) {
							temp +=
								this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) *
								cmR2[j()];
						}
						Kokkos::atomic_add(
							&values(_connEntries[2 + p][1 + k()][0][1 + i()]),
							this->_rate(gridIndex) * temp /
								(double)_productVolumes[p]);
					}

					// (d / dL_1^B)
					if (_reactantMomentIds[1][i()] != invalidIndex) {
						temp = this->_coefs(0, i() + 1, p + 2, k() + 1) * cR1;
						for (auto j : speciesRangeNoI) {
							temp +=
								this->_coefs(j() + 1, i() + 1, p + 2, k() + 1) *
								cmR1[j()];
						}
						Kokkos::atomic_add(
							&values(_connEntries[2 + p][1 + k()][1][1 + i()]),
							this->_rate(gridIndex) * temp /
								(double)_productVolumes[p]);
					}
				}
			}
		}
	}
}

KOKKOS_INLINE_FUNCTION
void
NEProductionReaction::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Initialize the concentrations that will be used in the loops
	auto cR1 = concentrations[_reactants[0]];
	Kokkos::Array<double, nMomentIds> cmR1;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[0][i()] == invalidIndex) {
			cmR1[i()] = 0.0;
		}
		else
			cmR1[i()] = concentrations[_reactantMomentIds[0][i()]];
	}
	auto cR2 = concentrations[_reactants[1]];
	Kokkos::Array<double, nMomentIds> cmR2;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] == invalidIndex) {
			cmR2[i()] = 0.0;
		}
		else
			cmR2[i()] = concentrations[_reactantMomentIds[1][i()]];
	}

	// Compute the partials for the 0th order moments
	// Compute the values (d / dL_0^A)
	double temp = this->_coefs(0, 0, 0, 0) * cR2;
	for (auto i : speciesRangeNoI) {
		temp += this->_coefs(0, i() + 1, 0, 0) * cmR2[i()];
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]),
		this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);
	// Second reactant
	if (_reactants[1] == _reactants[0])
		Kokkos::atomic_sub(&values(_connEntries[1][0][0][0]),
			this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex || prodId != _reactants[0]) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][0][0]),
			this->_rate(gridIndex) * temp / (double)_productVolumes[p]);
	}

	// Compute the values (d / dL_0^B)
	temp = this->_coefs(0, 0, 0, 0) * cR1;
	for (auto i : speciesRangeNoI) {
		temp += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()];
	}
	// First for the first reactant
	if (_reactants[1] == _reactants[0])
		Kokkos::atomic_sub(&values(_connEntries[0][0][1][0]),
			this->_rate(gridIndex) * temp / (double)_reactantVolumes[0]);
	// Second reactant
	Kokkos::atomic_sub(&values(_connEntries[1][0][1][0]),
		this->_rate(gridIndex) * temp / (double)_reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex || prodId != _reactants[1]) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][1][0]),
			this->_rate(gridIndex) * temp / (double)_productVolumes[p]);
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[0][k()] != invalidIndex) {
			// First for the first reactant
			for (auto i : speciesRangeNoI) {
				// (d / dL_1^A)
				temp = this->_coefs(i() + 1, 0, 0, k() + 1) * cR2;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(i() + 1, j() + 1, 0, k() + 1) * cmR2[j()];
				}
				if (k() == i())
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[0]);

				// (d / dL_1^B)
				temp = this->_coefs(0, i() + 1, 0, k() + 1) * cR1;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(j() + 1, i() + 1, 0, k() + 1) * cmR1[j()];
				}
				if (_reactantMomentIds[0][k()] == _reactantMomentIds[1][i()])
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][1][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[0]);
			}
		}

		if (_reactantMomentIds[1][k()] != invalidIndex) {
			// First for the second reactant
			for (auto i : speciesRangeNoI) {
				// (d / dL_1^A)
				temp = this->_coefs(i() + 1, 0, 1, k() + 1) * cR2;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(i() + 1, j() + 1, 1, k() + 1) * cmR2[j()];
				}
				if (_reactantMomentIds[1][k()] == _reactantMomentIds[0][i()])
					Kokkos::atomic_sub(
						&values(_connEntries[1][1 + k()][0][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[1]);

				// (d / dL_1^B)
				temp = this->_coefs(0, i() + 1, 1, k() + 1) * cR1;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(j() + 1, i() + 1, 1, k() + 1) * cmR1[j()];
				}
				if (k() == i())
					Kokkos::atomic_sub(
						&values(_connEntries[1][1 + k()][1][1 + i()]),
						this->_rate(gridIndex) * temp /
							(double)_reactantVolumes[1]);
			}
		}
	}

	// Loop on the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			if (_productMomentIds[p][k()] != invalidIndex) {
				for (auto i : speciesRangeNoI) {
					// (d / dL_1^A)
					temp = this->_coefs(i() + 1, 0, p + 2, k() + 1) * cR2;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) *
							cmR2[j()];
					}
					if (_productMomentIds[p][k()] == _reactantMomentIds[0][i()])
						Kokkos::atomic_add(
							&values(_connEntries[2 + p][1 + k()][0][1 + i()]),
							this->_rate(gridIndex) * temp /
								(double)_productVolumes[p]);

					// (d / dL_1^B)
					temp = this->_coefs(0, i() + 1, p + 2, k() + 1) * cR1;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, p + 2, k() + 1) *
							cmR1[j()];
					}
					if (_productMomentIds[p][k()] == _reactantMomentIds[1][i()])
						Kokkos::atomic_add(
							&values(_connEntries[2 + p][1 + k()][1][1 + i()]),
							this->_rate(gridIndex) * temp /
								(double)_productVolumes[p]);
				}
			}
		}
	}
}

KOKKOS_INLINE_FUNCTION
double
NEDissociationReaction::computeBindingEnergy()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);
	auto iFormation = this->_clusterData->getIFormationEnergy();
	auto vFormation = this->_clusterData->getVFormationEnergy();
	auto v2Formation = this->_clusterData->getV2FormationEnergy();
	auto xeFormation = this->_clusterData->getXeFormationEnergy();

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	Composition lo = clReg.getOrigin();
	Composition hi = clReg.getUpperLimitPoint();
	Composition prod1Comp = prod1Reg.getOrigin();
	Composition prod2Comp = prod2Reg.getOrigin();
	auto amtXe = (double)(lo[Species::Xe] + hi[Species::Xe] - 1) / 2.0;
	auto amtV = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;

	// Compute the delta_G = G(p1) + G(p2) - G(r)
	// Get the cluster ratio first
	double xeSD = amtXe / amtV;

	if (amtV == 0)
		return 5.0;

	// Fit parameters
	double fitParams[8] = {0.022098470861651, -0.328234592987279,
		1.9501944029492, -5.85771773178376, 9.18869701523602, -7.09118822290571,
		3.37190513184659, -1.03326651166341};

	if (prod1Comp.isOnAxis(Species::V) or prod2Comp.isOnAxis(Species::V)) {
		double xeSDPower = xeSD * xeSD;
		be = 0.0;
		for (auto i = 0; i < 8; i++) {
			be += (double)(i + 1) * fitParams[7 - i] * xeSDPower;
			xeSDPower *= xeSD;
		}

		if (prod1Comp.isOnAxis(Species::V) and prod1Comp[Species::V] == 2) {
			be *= 2.0;
			be += v2Formation;
		}
		else if (prod2Comp.isOnAxis(Species::V) and
			prod2Comp[Species::V] == 2) {
			be *= 2.0;
			be += v2Formation;
		}
		else {
			be += vFormation;
		}
	}
	if (prod1Comp.isOnAxis(Species::I) or prod2Comp.isOnAxis(Species::I)) {
		double xeSDPower = xeSD * xeSD;
		be = 0.0;
		for (auto i = 0; i < 8; i++) {
			be -= (double)(i + 1) * fitParams[7 - i] * xeSDPower;
			xeSDPower *= xeSD;
		}
		be += iFormation;
	}
	if (prod1Comp.isOnAxis(Species::Xe) or prod2Comp.isOnAxis(Species::Xe)) {
		double xeSDPower = xeSD;
		be = 0.0;
		for (auto i = 0; i < 8; i++) {
			be -= (double)(i + 2) * fitParams[7 - i] * xeSDPower;
			xeSDPower *= xeSD;
		}
		be += xeFormation;
	}

	return util::min(5.0, util::max(be, -1.0));
}

KOKKOS_INLINE_FUNCTION
double
NEDissociationReaction::getRateForProduction(IndexType gridIndex)
{
	return 0.0;
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return ne::getRate(cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

KOKKOS_INLINE_FUNCTION
double
NEDissociationReaction::computeRate(IndexType gridIndex, double time)
{
	double omega = this->_clusterData->atomicVolume();
	double T = this->_clusterData->temperature(gridIndex);
	constexpr double k_B = ::xolotl::core::kBoltzmann;

	// TODO: computeProductionRate should use products and not reactants
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	Composition cl0Comp = cl0.getRegion().getOrigin();
	Composition cl1Comp = cl1.getRegion().getOrigin();

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	//		if (cl0Comp[Species::Xe] > 0 and cl1Comp[Species::Xe] > 0)
	//			return 0.0;

	// Read the rates if available
	double rate = -1.0;
	if (this->_clusterData->extraData.fileClusterMap.exists(_products[0]) and
		this->_clusterData->extraData.fileClusterMap.exists(_products[1])) {
		auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
		rate = this->_clusterData->extraData.constantRates(
				   clusterMap.value_at(clusterMap.find(_products[0])),
				   clusterMap.value_at(clusterMap.find(_products[1])), 0) *
			exp(this->_deltaG0 / (k_B * T));
		if (this->_deltaG0 > 0.0) {
			rate = this->_clusterData->extraData.constantRates(
				clusterMap.value_at(clusterMap.find(_products[0])),
				clusterMap.value_at(clusterMap.find(_products[1])), 1);
		}
	}
	if (rate > 0.0) {
		return (1.0 / omega) * rate * 4.0 * ::xolotl::core::pi * r0 *
			::xolotl::core::zFactor;
	}

	double kPlus = getRateForProduction(gridIndex);
	double E_b = this->computeBindingEnergy();

	double kMinus = (1.0 / omega) * kPlus * exp(-E_b / (k_B * T));

	return kMinus;
}

KOKKOS_INLINE_FUNCTION
void
NEDissociationReaction::computeFlux(
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
	Kokkos::atomic_sub(&fluxes[_reactant], f / (double)_reactantVolume);
	Kokkos::atomic_add(&fluxes[_products[0]], f / (double)_productVolumes[0]);
	Kokkos::atomic_add(&fluxes[_products[1]], f / (double)_productVolumes[1]);

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
				&fluxes[_reactantMomentIds[k()]], f / (double)_reactantVolume);
		}

		// Now the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 1, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 1, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(&fluxes[_productMomentIds[0][k()]],
				f / (double)_productVolumes[0]);
		}

		// Finally the second product
		if (_productMomentIds[1][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 2, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 2, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(&fluxes[_productMomentIds[1][k()]],
				f / (double)_productVolumes[1]);
		}
	}
}

KOKKOS_INLINE_FUNCTION
void
NEDissociationReaction::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	using AmountType = typename NetworkType::AmountType;
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = this->_rate(gridIndex) / (double)_reactantVolume;
	// Compute the values
	Kokkos::atomic_sub(
		&values(_connEntries[0][0][0][0]), df * this->_coefs(0, 0, 0, 0));
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_sub(&values(_connEntries[0][0][0][1 + i()]),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the first product
	df = this->_rate(gridIndex) / (double)_productVolumes[0];
	Kokkos::atomic_add(
		&values(_connEntries[1][0][0][0]), df * this->_coefs(0, 0, 0, 0));

	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_add(&values(_connEntries[1][0][0][1 + i()]),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the second product
	df = this->_rate(gridIndex) / (double)_productVolumes[1];
	Kokkos::atomic_add(
		&values(_connEntries[2][0][0][0]), df * this->_coefs(0, 0, 0, 0));

	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_add(&values(_connEntries[2][0][0][1 + i()]),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = this->_rate(gridIndex) / (double)_reactantVolume;
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
		// For the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)_productVolumes[0];
			Kokkos::atomic_add(&values(_connEntries[1][1 + k()][0][0]),
				df * this->_coefs(0, 0, 1, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_add(
						&values(_connEntries[1][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 1, k() + 1));
				}
			}
		}
		// For the second product
		if (_productMomentIds[1][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)_productVolumes[1];
			Kokkos::atomic_add(&values(_connEntries[2][1 + k()][0][0]),
				df * this->_coefs(0, 0, 2, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_add(
						&values(_connEntries[2][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 2, k() + 1));
				}
			}
		}
	}
}

KOKKOS_INLINE_FUNCTION
void
NEDissociationReaction::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	using AmountType = typename NetworkType::AmountType;
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = this->_rate(gridIndex) / (double)_reactantVolume;
	// Compute the values
	Kokkos::atomic_sub(
		&values(_connEntries[0][0][0][0]), df * this->_coefs(0, 0, 0, 0));
	// For the first product
	df = this->_rate(gridIndex) / (double)_productVolumes[0];
	if (_products[0] == _reactant)
		Kokkos::atomic_add(
			&values(_connEntries[1][0][0][0]), df * this->_coefs(0, 0, 0, 0));

	// For the second product
	df = this->_rate(gridIndex) / (double)_productVolumes[1];
	if (_products[1] == _reactant)
		Kokkos::atomic_add(
			&values(_connEntries[2][0][0][0]), df * this->_coefs(0, 0, 0, 0));

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = this->_rate(gridIndex) / (double)_reactantVolume;
			// Compute the values
			for (auto i : speciesRangeNoI) {
				if (k() == i())
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 0, k() + 1));
			}
		}
		// For the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)_productVolumes[0];
			for (auto i : speciesRangeNoI) {
				if (_productMomentIds[0][k()] == _reactantMomentIds[i()])
					Kokkos::atomic_add(
						&values(_connEntries[1][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 1, k() + 1));
			}
		}
		// For the second product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)_productVolumes[1];
			for (auto i : speciesRangeNoI) {
				if (_productMomentIds[1][k()] == _reactantMomentIds[i()])
					Kokkos::atomic_add(
						&values(_connEntries[2][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 2, k() + 1));
			}
		}
	}
}

KOKKOS_INLINE_FUNCTION
double
NESinkReaction::getSinkBias()
{
	using Species = typename Superclass::Species;
	using Composition = typename Superclass::Composition;

	double bias = 1.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);

	auto clReg = cl.getRegion();
	if (clReg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		if (comp.isOnAxis(Species::I)) {
			bias = 1.1;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
NESinkReaction::getSinkStrength()
{
	// Not actually used
	//	return 1.0e-3;
	return 1.0e-4 * ::xolotl::core::zFactor;
}

KOKKOS_INLINE_FUNCTION
void
NESinkReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
	double omega = this->_clusterData->atomicVolume();
	auto rate = ::xolotl::core::zFactor * getSinkBias() *
		this->_clusterData->extraData.constantRates(
			clusterMap.value_at(clusterMap.find(_reactant)),
			this->_clusterData->extraData.constantRates.extent(0), 0);
	// delta G_0 is assumed to be always negative
	Kokkos::atomic_add(&fluxes(_reactant),
		rate *
			(exp(this->_deltaG0 /
				 (::xolotl::core::kBoltzmann *
					 this->_clusterData->temperature(gridIndex))) /
					omega -
				concentrations(_reactant)));

	return;
}

KOKKOS_INLINE_FUNCTION
void
NESinkReaction::computePartialDerivatives(ConcentrationsView concentrations,
	Kokkos::View<double*> values, IndexType gridIndex)
{
	auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
	auto rate = this->_clusterData->extraData.constantRates(
					clusterMap.value_at(clusterMap.find(_reactant)),
					this->_clusterData->extraData.constantRates.extent(0), 0) *
		::xolotl::core::zFactor * getSinkBias();
	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), rate);
	return;
}

KOKKOS_INLINE_FUNCTION
void
NESinkReaction::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	auto& clusterMap = this->_clusterData->extraData.fileClusterMap;
	auto rate = this->_clusterData->extraData.constantRates(
					clusterMap.value_at(clusterMap.find(_reactant)),
					this->_clusterData->extraData.constantRates.extent(0), 0) *
		::xolotl::core::zFactor * getSinkBias();
	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), rate);
	return;
}
} // namespace network
} // namespace core
} // namespace xolotl
