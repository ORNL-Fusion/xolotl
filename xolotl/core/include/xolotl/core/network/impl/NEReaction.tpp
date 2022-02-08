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
KOKKOS_INLINE_FUNCTION
double
NEProductionReaction::computeRate(IndexType gridIndex)
{
	// Read the rates if available
	auto rate = this->_clusterData->extraData.constantRates(
		_reactants[0], _reactants[1]);
	if (rate > 0) {
		return rate;
	}

	auto cl0 = this->_clusterData->getCluster(_reactants[0]);
	auto cl1 = this->_clusterData->getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	rate = getRateForProduction(
		cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);

	//	constexpr auto speciesRange = NetworkType::getSpeciesRange();
	//	auto cl0Reg = cl0.getRegion(), cl1Reg = cl1.getRegion(),
	//		 prodReg = this->_clusterData.getCluster(_products[0]).getRegion();
	//	for (auto i : speciesRange) {
	//		std::cout << cl0Reg[i()].begin() << " ";
	//	}
	//	std::cout << "+ ";
	//	for (auto i : speciesRange) {
	//		std::cout << cl1Reg[i()].begin() << " ";
	//	}
	//	std::cout << "-> ";
	//	for (auto i : speciesRange) {
	//		std::cout << prodReg[i()].begin() << " ";
	//	}
	//	std::cout << std::endl;

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
		// Compute the flux for the 0th order moments
		double f = ::xolotl::core::uConcentration;
		f *= this->_clusterData->extraData.constantRates(
				 _reactants[0], _reactants[1]) *
			std::exp(this->_deltaG0 / ::xolotl::core::kBoltzmann *
				this->_clusterData->temperature(gridIndex));

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

	constexpr double xeTrapTable[9] = {
		0.0, 4.31, 2.90, 2.02, 1.09, 0.58, 0.13, -0.25, -0.59};

	double be = 5.0;

	auto cl = this->_clusterData->getCluster(this->_reactant);
	auto prod1 = this->_clusterData->getCluster(this->_products[0]);
	auto prod2 = this->_clusterData->getCluster(this->_products[1]);

	auto clReg = cl.getRegion();
	auto prod1Reg = prod1.getRegion();
	auto prod2Reg = prod2.getRegion();
	if (clReg.isSimplex() && prod1Reg.isSimplex() && prod2Reg.isSimplex()) {
		Composition comp = clReg.getOrigin();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		if (comp.isOnAxis(Species::Xe)) {
			if (prod1Comp.isOnAxis(Species::Xe) ||
				prod2Comp.isOnAxis(Species::Xe)) {
				if (comp[Species::Xe] == 2)
					be = 0.5;
				else
					be = 1.0;
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				be = xeTrapTable[comp[Species::Xe]];
			}
		}
		else if (comp.isOnAxis(Species::V)) {
			auto size = comp[Species::V];
			be = 1.73 -
				2.59 *
					(pow((double)size, 2.0 / 3.0) -
						pow((double)size - 1.0, 2.0 / 3.0));
		}
		else if (comp.isOnAxis(Species::I)) {
			// Nothing
		}
		else {
			// XeV
			if (prod1Comp.isOnAxis(Species::V) ||
				prod2Comp.isOnAxis(Species::V)) {
				auto amtXe = comp[Species::Xe], amtV = comp[Species::V];
				be = 1.73 -
					2.59 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) +
					2.5 * log(1.0 + ((double)amtXe / (double)amtV));
			}
			if (prod1Comp.isOnAxis(Species::I) ||
				prod2Comp.isOnAxis(Species::I)) {
				auto amtXe = comp[Species::Xe], amtV = comp[Species::V];
				be = 4.88 +
					2.59 *
						(pow((double)amtV, 2.0 / 3.0) -
							pow((double)amtV - 1.0, 2.0 / 3.0)) -
					2.5 * log(1.0 + ((double)amtXe / (double)amtV));
			}
		}
	}
	else {
		Composition lo = clReg.getOrigin();
		Composition hi = clReg.getUpperLimitPoint();
		Composition prod1Comp = prod1Reg.getOrigin();
		Composition prod2Comp = prod2Reg.getOrigin();
		// XeV
		auto amtXe = (double)(lo[Species::Xe] + hi[Species::Xe] - 1) / 2.0;
		auto amtV = (double)(lo[Species::V] + hi[Species::V] - 1) / 2.0;
		if (prod1Comp.isOnAxis(Species::V) || prod2Comp.isOnAxis(Species::V)) {
			be = 1.73 -
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) +
				2.5 * log(1.0 + (amtXe / amtV));
		}
		if (prod1Comp.isOnAxis(Species::I) || prod2Comp.isOnAxis(Species::I)) {
			be = 4.88 +
				2.59 * (pow(amtV, 2.0 / 3.0) - pow(amtV - 1.0, 2.0 / 3.0)) -
				2.5 * log(1.0 + (amtXe / amtV));
		}
	}

	return util::min(5.0, util::max(be, -5.0));
}

KOKKOS_INLINE_FUNCTION
double
NEDissociationReaction::computeRate(IndexType gridIndex)
{
	double omega = this->_clusterData->atomicVolume();
	double T = this->_clusterData->temperature(gridIndex);
	constexpr double k_B = ::xolotl::core::kBoltzmann;

	// Read the rates if available
	auto rate =
		this->_clusterData->extraData.constantRates(_products[0], _products[1]);
	if (rate > 0) {
		return rate * std::exp(this->_deltaG0 / (k_B * T));
	}

	// TODO: computeProductionRate should use products and not reactants
	auto cl0 = this->_clusterData->getCluster(_products[0]);
	auto cl1 = this->_clusterData->getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	double kPlus = getRateForProduction(
		cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
	double E_b = this->asDerived()->computeBindingEnergy();

	//	if (E_b < 1.0) {
	//	auto cl0Reg = cl0.getRegion().getOrigin(), cl1Reg =
	// cl1.getRegion().getOrigin(), prod0Reg =
	// this->_clusterData.getCluster(_reactant).getRegion().getOrigin();
	//
	//	constexpr auto speciesRange = NetworkType::getSpeciesRange();
	//	for (auto j : speciesRange) {
	//				std::cout << cl0Reg[j()] << " ";
	//			}
	//	std::cout << r0 << std::endl << " + " << std::endl;
	//	for (auto j : speciesRange) {
	//				std::cout << cl1Reg[j()] << " ";
	//			}
	//	std::cout << r1 << std::endl << " -> " << std::endl;
	//	for (auto j : speciesRange) {
	//				std::cout << prod0Reg[j()] << " ";
	//			}
	//	std::cout << E_b << std::endl << "." << std::endl;
	//	}

	double kMinus = (1.0 / omega) * kPlus * std::exp(-E_b / (k_B * T));

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
			bias = 1.05;
		}
	}

	return bias;
}

KOKKOS_INLINE_FUNCTION
double
NESinkReaction::getSinkStrength()
{
	auto cl = this->_clusterData->getCluster(this->_reactant);
	double r = cl.getReactionRadius();
	double latticeParameter = this->_clusterData->latticeParameter();
	double r0 = latticeParameter * 0.5 * sqrt(2.0);
	double rho = 0.0003;
	constexpr double pi = ::xolotl::core::pi;

	double strength = -4.0 * pi * rho / log(pi * rho * (r + r0) * (r + r0));

	return strength;
}

KOKKOS_INLINE_FUNCTION
void
NESinkReaction::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	Kokkos::atomic_add(&fluxes(_reactant),
		this->_clusterData->extraData.constantRates(
			_reactant, this->_clusterData->numClusters) *
				::xolotl::core::uConcentration *
				std::exp(this->_deltaG0 /
					(::xolotl::core::kBoltzmann *
						this->_clusterData->temperature(gridIndex))) -
			this->_clusterData->extraData.constantRates(
				_reactant, this->_clusterData->numClusters) *
				concentrations(_reactant));
}

KOKKOS_INLINE_FUNCTION
void
NESinkReaction::computePartialDerivatives(ConcentrationsView concentrations,
	Kokkos::View<double*> values, IndexType gridIndex)
{
	// Read the rates if available
	auto rate = this->_clusterData->extraData.constantRates(
		_reactant, this->_clusterData->numClusters);

	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), rate);
}

KOKKOS_INLINE_FUNCTION
void
NESinkReaction::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	IndexType gridIndex)
{
	// Read the rates if available
	auto rate = this->_clusterData->extraData.constantRates(
		_reactant, this->_clusterData->numClusters);

	Kokkos::atomic_sub(&values(_connEntries[0][0][0][0]), rate);
}
} // namespace network
} // namespace core
} // namespace xolotl
