#pragma once

#include <algorithm>
#include <array>

#include <plsm/EnumIndexed.h>

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/ReactionUtility.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION
Reaction<TNetwork, TDerived>::Reaction(ReactionDataRef reactionData,
	const ClusterData& clusterData, IndexType reactionId) :
	_clusterData(&clusterData),
	_reactionId(reactionId),
	_rate(reactionData.getRates(reactionId)),
	_widths(reactionData.getWidths(reactionId)),
	_coefs(reactionData.getCoefficients(reactionId)),
	_deltaG0(0.0)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::updateData(
	ReactionDataRef reactionData, const ClusterData& clusterData)
{
	_clusterData = &clusterData;
	_rate = reactionData.getRates(_reactionId);
	_constantRates = reactionData.getConstantRates(_reactionId);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::getRateEntries(ReactionDataRef reactionData)
{
	_rateEntries = reactionData.getRateEntries(_reactionId);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
Reaction<TNetwork, TDerived>::computeOverlap(const ReflectedRegion& cl1RR,
	const ReflectedRegion& cl2RR, const ReflectedRegion& pr1RR,
	const ReflectedRegion& pr2RR)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	double nOverlap = 1.0;
	for (auto i : speciesRangeNoI) {
		// The width is the subset of the tiles for which the
		// reaction is possible
		// For instance if we have X_1 + X_[3,5) ⇄ X_[5,7)
		// It is only possible for 4 within X_[3,5) and 5 within X_[5,7)
		// so the width is 1
		// More complicated with X_[3,5) + X_[5,7) ⇄ X_[9,11)
		// 3+6, 4+5, 4+6, width is 3

		for (auto m : makeIntervalRange(pr2RR[i()]))
			for (auto l : makeIntervalRange(cl2RR[i()])) {
				_widths(i()) += util::max(0.0,
					util::min(
						cl1RR[i()].end() - 1 + l, pr1RR[i()].end() - 1 + m) -
						util::max(
							cl1RR[i()].begin() + l, pr1RR[i()].begin() + m) +
						1.0);
			}

		nOverlap *= _widths(i());
	}

	//	if (nOverlap <= 0) {
	//	std::cout << "first reactant: ";
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << cl1RR[i()].begin() << ", ";
	//	}
	//	std::cout << std::endl;
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << cl1RR[i()].end() - 1 << ", ";
	//	}
	//	std::cout << std::endl << "second reactant: ";
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << cl2RR[i()].begin() << ", ";
	//	}
	//	std::cout << std::endl;
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << cl2RR[i()].end() - 1 << ", ";
	//	}
	//	std::cout << std::endl << "product: ";
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << pr1RR[i()].begin() << ", ";
	//	}
	//	std::cout << std::endl;
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << pr1RR[i()].end() - 1 << ", ";
	//	}
	//	std::cout << std::endl << "second product: ";
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << pr2RR[i()].begin() << ", ";
	//	}
	//	std::cout << std::endl;
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << pr2RR[i()].end() - 1 << ", ";
	//	}
	//	std::cout << std::endl;
	//	std::cout << "Overlap: " << nOverlap << std::endl;
	//	std::cout << "Widths: ";
	//	for (auto i : speciesRangeNoI) {
	//		std::cout << _widths(i()) << ", ";
	//	}
	//	std::cout << std::endl;
	//	}

	assert(nOverlap > 0);

	return nOverlap;
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION
ProductionReaction<TNetwork, TDerived>::ProductionReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1,
	IndexType cluster2, IndexType cluster3) :
	Superclass(reactionData, clusterData, reactionId),
	_reactants({cluster0, cluster1}),
	_products({cluster2, cluster3})
{
	for (auto i : {0, 1}) {
		this->copyMomentIds(_reactants[i], _reactantMomentIds[i]);
		this->copyMomentIds(_products[i], _productMomentIds[i]);
	}

	// static
	const auto dummyRegion = Region(Composition{});

	const auto& cl1Reg =
		this->_clusterData->getCluster(_reactants[0]).getRegion();
	const auto& cl2Reg =
		this->_clusterData->getCluster(_reactants[1]).getRegion();
	const auto& pr1Reg = (_products[0] == invalidIndex) ?
		dummyRegion :
		this->_clusterData->getCluster(_products[0]).getRegion();
	const auto& pr2Reg = (_products[1] == invalidIndex) ?
		dummyRegion :
		this->_clusterData->getCluster(_products[1]).getRegion();

	_reactantVolumes = {cl1Reg.volume(), cl2Reg.volume()};
	_productVolumes = {pr1Reg.volume(), pr2Reg.volume()};

	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION
ProductionReaction<TNetwork, TDerived>::ProductionReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	ProductionReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2,
		clusterSet.cluster3)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeCoefficients()
{
	// static
	const auto dummyRegion = Region(Composition{});

	// Find the overlap for this reaction
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	const auto& cl1Reg =
		this->_clusterData->getCluster(_reactants[0]).getRegion();
	const auto& cl2Reg =
		this->_clusterData->getCluster(_reactants[1]).getRegion();
	const auto& pr1Reg = (_products[0] == invalidIndex) ?
		dummyRegion :
		this->_clusterData->getCluster(_products[0]).getRegion();
	const auto& pr2Reg = (_products[1] == invalidIndex) ?
		dummyRegion :
		this->_clusterData->getCluster(_products[1]).getRegion();
	const auto& cl1Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			cl1Reg);
	const auto& cl2Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			cl2Reg);

	// Initialize the reflected regions
	auto rRegions = detail::updateReflectedRegionsForCoefs<nMomentIds>(
		cl1Reg, cl2Reg, pr1Reg, pr2Reg);
	auto cl1RR = rRegions[0];
	auto cl2RR = rRegions[1];
	auto pr1RR = rRegions[2];
	auto pr2RR = rRegions[3];

	// If there is no product the overlap is 1
	double nOverlap = 1.0;
	// No product case
	if (_products[0] == invalidIndex && _products[1] == invalidIndex) {
		for (auto i : speciesRangeNoI) {
			this->_widths[i()] = 1.0;
		}
	}
	// General case
	else
		nOverlap = this->computeOverlap(cl1RR, cl2RR, pr1RR, pr2RR);

	this->_coefs(0, 0, 0, 0) = nOverlap;
	for (auto i : speciesRangeNoI) {
		// First order sum on the first reactant
		auto factor = nOverlap / this->_widths[i()];
		this->_coefs(i() + 1, 0, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), cl1RR, cl2RR, pr1RR, pr2RR);

		this->_coefs(0, 0, 0, i() + 1) =
			this->_coefs(i() + 1, 0, 0, 0) / cl1Disp[i()];

		// First order sum on the second reactant
		this->_coefs(0, i() + 1, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), cl2RR, cl1RR, pr1RR, pr2RR);

		this->_coefs(0, 0, 1, i() + 1) =
			this->_coefs(0, i() + 1, 0, 0) / cl2Disp[i()];

		// Loop on the potential products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}

			// Get the regions in the right order
			const auto& thisRR = (prodId == _products[0]) ? pr1RR : pr2RR;
			const auto& otherRR = (prodId == _products[0]) ? pr2RR : pr1RR;
			// Get the dispersion
			const auto& thisDispersion = (prodId == _products[0]) ?
				detail::getReflectedDispersionForCoefs<
					NetworkType::Traits::numSpecies>(pr1Reg) :
				detail::getReflectedDispersionForCoefs<
					NetworkType::Traits::numSpecies>(pr2Reg);

			// First order sum on the other product (p+2) because 0 and 1 are
			// used for reactants)
			this->_coefs(0, 0, p + 2, i() + 1) = factor *
				detail::computeFirstOrderSum(
					i(), thisRR, otherRR, cl2RR, cl1RR) /
				thisDispersion[i()];

			// Products first moments
			for (auto k : speciesRangeNoI) {
				// Second order sum
				if (k == i) {
					this->_coefs(i() + 1, 0, p + 2, k() + 1) = factor *
						detail::computeSecondOrderOffsetSum(
							i(), cl1RR, cl2RR, thisRR, otherRR) /
						thisDispersion[i()];

					this->_coefs(0, i() + 1, p + 2, k() + 1) = factor *
						detail::computeSecondOrderOffsetSum(
							i(), cl2RR, cl1RR, thisRR, otherRR) /
						thisDispersion[i()];
				}
				else {
					this->_coefs(i() + 1, 0, p + 2, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, 0, p + 2, k() + 1) / nOverlap;

					this->_coefs(0, i() + 1, p + 2, k() + 1) =
						this->_coefs(0, i() + 1, 0, 0) *
						this->_coefs(0, 0, p + 2, k() + 1) / nOverlap;
				}
			}
		}
	}

	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		for (auto j : speciesRangeNoI) {
			// Second order sum
			if (i == j) {
				for (double m : makeIntervalRange(pr2RR[j()]))
					for (double l : makeIntervalRange(cl1RR[j()])) {
						this->_coefs(i() + 1, j() + 1, 0, 0) +=
							(l -
								static_cast<double>(
									cl1RR[j()].end() - 1 + cl1RR[j()].begin()) /
									2.0) *
							factor *
							util::firstOrderSum(
								util::max(pr1RR[j()].begin() + m - l,
									static_cast<double>(cl2RR[j()].begin())),
								util::min(pr1RR[j()].end() - 1 + m - l,
									static_cast<double>(cl2RR[j()].end() - 1)),
								static_cast<double>(
									cl2RR[j()].end() - 1 + cl2RR[j()].begin()) /
									2.0);
					}
			}
			else {
				this->_coefs(i() + 1, j() + 1, 0, 0) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(0, j() + 1, 0, 0) / nOverlap;
			}
		}

		factor = nOverlap / this->_widths[i()];

		// First reactant first moments
		for (auto k : speciesRangeNoI) {
			if (k == i) {
				this->_coefs(i() + 1, 0, 0, k() + 1) = factor *
					detail::computeSecondOrderSum(
						i(), cl1RR, cl2RR, pr1RR, pr2RR) /
					cl1Disp[i()];
				this->_coefs(0, i() + 1, 1, k() + 1) = factor *
					detail::computeSecondOrderSum(
						i(), cl2RR, cl1RR, pr1RR, pr2RR) /
					cl2Disp[i()];
			}
			else {
				this->_coefs(i() + 1, 0, 0, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(k() + 1, 0, 0, 0) / (nOverlap * cl1Disp[k()]);
				this->_coefs(0, i() + 1, 1, k() + 1) =
					this->_coefs(0, i() + 1, 0, 0) *
					this->_coefs(0, k() + 1, 0, 0) / (nOverlap * cl2Disp[k()]);
			}

			this->_coefs(0, i() + 1, 0, k() + 1) =
				this->_coefs(k() + 1, i() + 1, 0, 0) / cl1Disp[k()];

			// Second reactant partial derivatives
			this->_coefs(i() + 1, 0, 1, k() + 1) =
				this->_coefs(i() + 1, k() + 1, 0, 0) / cl2Disp[k()];
		}

		// Now we loop over the 2 dimensions of the coefs to compute all
		// the possible sums over distances for the flux
		factor = nOverlap / this->_widths[i()];
		for (auto j : speciesRangeNoI) {
			// Now we deal with the coefficients needed for the
			// first moments
			// Let's start with the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}

				// Get the regions in the right order
				const auto& thisRR = (prodId == _products[0]) ? pr1RR : pr2RR;
				const auto& otherRR = (prodId == _products[0]) ? pr2RR : pr1RR;
				// Get the dispersion
				const auto& thisDispersion = (prodId == _products[0]) ?
					detail::getReflectedDispersionForCoefs<
						NetworkType::Traits::numSpecies>(pr1Reg) :
					detail::getReflectedDispersionForCoefs<
						NetworkType::Traits::numSpecies>(pr2Reg);

				for (auto k : speciesRangeNoI) {
					// Third order sum
					if (i == j && j == k) {
						for (double m : makeIntervalRange(otherRR[i()]))
							for (double l : makeIntervalRange(cl1RR[i()])) {
								this->_coefs(
									i() + 1, j() + 1, p + 2, k() + 1) +=
									(l -
										static_cast<double>(cl1RR[i()].end() -
											1 + cl1RR[i()].begin()) /
											2.0) *
									factor *
									util::secondOrderOffsetSum(
										util::max(thisRR[i()].begin() + m - l,
											static_cast<double>(
												cl2RR[i()].begin())),
										util::min(thisRR[i()].end() - 1 + m - l,
											static_cast<double>(
												cl2RR[i()].end() - 1)),
										static_cast<double>(cl2RR[i()].end() -
											1 + cl2RR[i()].begin()) /
											2.0,
										static_cast<double>(thisRR[i()].end() -
											1 + thisRR[i()].begin()) /
											2.0,
										l - m);
							}
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) /=
							thisDispersion[k()];
					}
					else if (j == k) {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(i() + 1, 0, 0, 0) *
							this->_coefs(0, j() + 1, p + 2, k() + 1) / nOverlap;
					}
					else if (i == k) {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(0, j() + 1, 0, 0) *
							this->_coefs(i() + 1, 0, p + 2, k() + 1) / nOverlap;
					}
					else if (i == j) {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(0, 0, p + 2, k() + 1) *
							this->_coefs(i() + 1, j() + 1, 0, 0) / nOverlap;
					}
					else {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(i() + 1, 0, 0, 0) *
							this->_coefs(0, j() + 1, 0, 0) *
							this->_coefs(0, 0, p + 2, k() + 1) /
							(nOverlap * nOverlap);
					}
				}
			}

			// Let's take care of the first reactant first moments
			for (auto k : speciesRangeNoI) {
				// Third order sum
				if (i == j && j == k) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) = factor *
						detail::computeThirdOrderSum(
							i(), cl2RR, cl1RR, pr1RR, pr2RR) /
						cl1Disp[i()];
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) = factor *
						detail::computeThirdOrderSum(
							i(), cl1RR, cl2RR, pr1RR, pr2RR) /
						cl2Disp[i()];
				}
				else if (i == k) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(i() + 1, 0, 0, k() + 1) / nOverlap;
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(i() + 1, 0, 1, k() + 1) / nOverlap;
				}
				else if (j == k) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 0, k() + 1) / nOverlap;
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 1, k() + 1) / nOverlap;
				}
				else if (i == j) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(0, 0, 0, k() + 1) *
						this->_coefs(i() + 1, j() + 1, 0, 0) / nOverlap;
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(0, 0, 1, k() + 1) *
						this->_coefs(i() + 1, j() + 1, 0, 0) / nOverlap;
				}
				else {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(k() + 1, 0, 0, 0) /
						(nOverlap * nOverlap * cl1Disp[k()]);
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(0, k() + 1, 0, 0) /
						(nOverlap * nOverlap * cl2Disp[k()]);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
ProductionReaction<TNetwork, TDerived>::computeRate(
	IndexType gridIndex, double time)
{
	return this->asDerived()->getRateForProduction(gridIndex);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	// Each reactant connects with all the reactants
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

	// Reactant 2 with reactant 1
	this->addConnectivity(_reactants[1], _reactants[0], connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[0][i()] != invalidIndex) {
			this->addConnectivity(
				_reactants[1], _reactantMomentIds[0][i()], connectivity);
		}
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			this->addConnectivity(
				_reactantMomentIds[1][i()], _reactants[0], connectivity);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[0][j()] != invalidIndex) {
					this->addConnectivity(_reactantMomentIds[1][i()],
						_reactantMomentIds[0][j()], connectivity);
				}
			}
		}
	}
	// Reactant 1 with reactant 2
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
	// Reactant 2 with reactant 2
	this->addConnectivity(_reactants[1], _reactants[1], connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			this->addConnectivity(
				_reactants[1], _reactantMomentIds[1][i()], connectivity);
			this->addConnectivity(
				_reactantMomentIds[1][i()], _reactants[1], connectivity);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[1][j()] != invalidIndex) {
					this->addConnectivity(_reactantMomentIds[1][i()],
						_reactantMomentIds[1][j()], connectivity);
				}
			}
		}
	}
	// Each product connects with all the reactants
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		auto prod = this->_clusterData->getCluster(prodId);
		const auto& prodReg = prod.getRegion();

		// With reactant 1
		this->addConnectivity(prodId, _reactants[0], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				this->addConnectivity(
					prodId, _reactantMomentIds[0][i()], connectivity);
			}
			if (_productMomentIds[p][i()] != invalidIndex) {
				this->addConnectivity(
					_productMomentIds[p][i()], _reactants[0], connectivity);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[0][j()] != invalidIndex) {
						this->addConnectivity(_productMomentIds[p][i()],
							_reactantMomentIds[0][j()], connectivity);
					}
				}
			}
		}
		// With reactant 2
		this->addConnectivity(prodId, _reactants[1], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] != invalidIndex) {
				this->addConnectivity(
					prodId, _reactantMomentIds[1][i()], connectivity);
			}
			if (_productMomentIds[p][i()] != invalidIndex) {
				this->addConnectivity(
					_productMomentIds[p][i()], _reactants[1], connectivity);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[1][j()] != invalidIndex) {
						this->addConnectivity(_productMomentIds[p][i()],
							_reactantMomentIds[1][j()], connectivity);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
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
	// Reactant 2 with reactant 2
	this->addConnectivity(_reactants[1], _reactants[1], connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (i() == j())
					this->addConnectivity(_reactantMomentIds[1][i()],
						_reactantMomentIds[1][j()], connectivity);
			}
		}
	}
	// Each product connects with all the reactants
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		auto prod = this->_clusterData->getCluster(prodId);
		const auto& prodReg = prod.getRegion();

		// With reactant 1
		if (prodId == _reactants[0])
			this->addConnectivity(prodId, _reactants[0], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_productMomentIds[p][i()] != invalidIndex) {
				for (auto j : speciesRangeNoI) {
					if (_productMomentIds[p][i()] == _reactantMomentIds[0][j()])
						this->addConnectivity(_productMomentIds[p][i()],
							_reactantMomentIds[0][j()], connectivity);
				}
			}
		}
		// With reactant 2
		if (prodId == _reactants[1])
			this->addConnectivity(prodId, _reactants[1], connectivity);
		for (auto i : speciesRangeNoI) {
			if (_productMomentIds[p][i()] != invalidIndex) {
				for (auto j : speciesRangeNoI) {
					if (_productMomentIds[p][i()] == _reactantMomentIds[1][j()])
						this->addConnectivity(_productMomentIds[p][i()],
							_reactantMomentIds[1][j()], connectivity);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
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

	Kokkos::atomic_sub(&fluxes[_reactants[0]], f / _reactantVolumes[0]);
	Kokkos::atomic_sub(&fluxes[_reactants[1]], f / _reactantVolumes[1]);

	IndexType p = 0;
	for (auto prodId : _products) {
		if (prodId == invalidIndex) {
			continue;
		}
		Kokkos::atomic_add(&fluxes[prodId], f / _productVolumes[p]);
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
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[0][k()]], f / _reactantVolumes[0]);
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
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[1][k()]], f / _reactantVolumes[1]);
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
				Kokkos::atomic_add(
					&fluxes[_productMomentIds[p][k()]], f / _productVolumes[p]);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computePartialDerivatives(
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
		this->_rate(gridIndex) * temp / _reactantVolumes[0]);
	// Second reactant
	Kokkos::atomic_sub(&values(_connEntries[1][0][0][0]),
		this->_rate(gridIndex) * temp / _reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][0][0]),
			this->_rate(gridIndex) * temp / _productVolumes[p]);
	}

	// Compute the values (d / dL_0^B)
	temp = this->_coefs(0, 0, 0, 0) * cR1;
	for (auto i : speciesRangeNoI) {
		temp += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()];
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(_connEntries[0][0][1][0]),
		this->_rate(gridIndex) * temp / _reactantVolumes[0]);
	// Second reactant
	Kokkos::atomic_sub(&values(_connEntries[1][0][1][0]),
		this->_rate(gridIndex) * temp / _reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][1][0]),
			this->_rate(gridIndex) * temp / _productVolumes[p]);
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
				this->_rate(gridIndex) * temp / _reactantVolumes[0]);
			// second reactant
			Kokkos::atomic_sub(&values(_connEntries[1][0][0][1 + i()]),
				this->_rate(gridIndex) * temp / _reactantVolumes[1]);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][0][0][1 + i()]),
					this->_rate(gridIndex) * temp / _productVolumes[p]);
			}
		}

		// (d / dL_1^B)
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			temp = this->_coefs(0, i() + 1, 0, 0) * cR1;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, i() + 1, 0, 0) * cmR1[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[0][0][1][1 + i()]),
				this->_rate(gridIndex) * temp / _reactantVolumes[0]);
			Kokkos::atomic_sub(&values(_connEntries[1][0][1][1 + i()]),
				this->_rate(gridIndex) * temp / _reactantVolumes[1]);
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][0][1][1 + i()]),
					this->_rate(gridIndex) * temp / _productVolumes[p]);
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
				this->_rate(gridIndex) * temp / _reactantVolumes[0]);

			// (d / dL_0^B)
			temp = this->_coefs(0, 0, 0, k() + 1) * cR1;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, 0, 0, k() + 1) * cmR1[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[0][1 + k()][1][0]),
				this->_rate(gridIndex) * temp / _reactantVolumes[0]);

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
						this->_rate(gridIndex) * temp / _reactantVolumes[0]);
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
						this->_rate(gridIndex) * temp / _reactantVolumes[0]);
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
				this->_rate(gridIndex) * temp / _reactantVolumes[1]);

			// (d / dL_0^B)
			temp = this->_coefs(0, 0, 1, k() + 1) * cR1;
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, 0, 1, k() + 1) * cmR1[j()];
			}
			Kokkos::atomic_sub(&values(_connEntries[1][1 + k()][1][0]),
				this->_rate(gridIndex) * temp / _reactantVolumes[1]);

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
						this->_rate(gridIndex) * temp / _reactantVolumes[1]);
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
						this->_rate(gridIndex) * temp / _reactantVolumes[1]);
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
					this->_rate(gridIndex) * temp / _productVolumes[p]);

				// (d / dL_0^B)
				temp = this->_coefs(0, 0, p + 2, k() + 1) * cR1;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(j() + 1, 0, p + 2, k() + 1) * cmR1[j()];
				}
				Kokkos::atomic_add(&values(_connEntries[2 + p][1 + k()][1][0]),
					this->_rate(gridIndex) * temp / _productVolumes[p]);

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
							this->_rate(gridIndex) * temp / _productVolumes[p]);
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
							this->_rate(gridIndex) * temp / _productVolumes[p]);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
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
		this->_rate(gridIndex) * temp / _reactantVolumes[0]);
	// Second reactant
	if (_reactants[1] == _reactants[0])
		Kokkos::atomic_sub(&values(_connEntries[1][0][0][0]),
			this->_rate(gridIndex) * temp / _reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex || prodId != _reactants[0]) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][0][0]),
			this->_rate(gridIndex) * temp / _productVolumes[p]);
	}

	// Compute the values (d / dL_0^B)
	temp = this->_coefs(0, 0, 0, 0) * cR1;
	for (auto i : speciesRangeNoI) {
		temp += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()];
	}
	// First for the first reactant
	if (_reactants[1] == _reactants[0])
		Kokkos::atomic_sub(&values(_connEntries[0][0][1][0]),
			this->_rate(gridIndex) * temp / _reactantVolumes[0]);
	// Second reactant
	Kokkos::atomic_sub(&values(_connEntries[1][0][1][0]),
		this->_rate(gridIndex) * temp / _reactantVolumes[1]);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex || prodId != _reactants[1]) {
			continue;
		}
		Kokkos::atomic_add(&values(_connEntries[2 + p][0][1][0]),
			this->_rate(gridIndex) * temp / _productVolumes[p]);
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
						this->_rate(gridIndex) * temp / _reactantVolumes[0]);

				// (d / dL_1^B)
				temp = this->_coefs(0, i() + 1, 0, k() + 1) * cR1;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(j() + 1, i() + 1, 0, k() + 1) * cmR1[j()];
				}
				if (_reactantMomentIds[0][k()] == _reactantMomentIds[1][i()])
					Kokkos::atomic_sub(
						&values(_connEntries[0][1 + k()][1][1 + i()]),
						this->_rate(gridIndex) * temp / _reactantVolumes[0]);
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
						this->_rate(gridIndex) * temp / _reactantVolumes[1]);

				// (d / dL_1^B)
				temp = this->_coefs(0, i() + 1, 1, k() + 1) * cR1;
				for (auto j : speciesRangeNoI) {
					temp +=
						this->_coefs(j() + 1, i() + 1, 1, k() + 1) * cmR1[j()];
				}
				if (k() == i())
					Kokkos::atomic_sub(
						&values(_connEntries[1][1 + k()][1][1 + i()]),
						this->_rate(gridIndex) * temp / _reactantVolumes[1]);
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
							this->_rate(gridIndex) * temp / _productVolumes[p]);

					// (d / dL_1^B)
					temp = this->_coefs(0, i() + 1, p + 2, k() + 1) * cR1;
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, p + 2, k() + 1) *
							cmR1[j()];
					}
					if (_productMomentIds[p][k()] == _reactantMomentIds[1][i()])
						Kokkos::atomic_add(
							&values(_connEntries[2 + p][1 + k()][1][1 + i()]),
							this->_rate(gridIndex) * temp / _productVolumes[p]);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeConstantRates(
	ConcentrationsView concentrations, RatesView rates, BelongingView isInSub,
	IndexType subId, IndexType gridIndex)
{
	// Check products
	bool productInSub = false;
	AmountType nProd = 0;
	for (auto prodId : _products) {
		if (prodId == invalidIndex) {
			continue;
		}
		nProd++;
		if (isInSub[prodId])
			productInSub = true;
	}
	// Only consider specific cases
	if (not isInSub[_reactants[0]] and not isInSub[_reactants[1]]) {
		if (nProd == 0)
			return;
		if (nProd > 0 && not productInSub)
			return;
	}
	if (isInSub[_reactants[0]] and isInSub[_reactants[1]]) {
		if (nProd == 0)
			return;
		if (nProd > 0 && productInSub)
			return;
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

	auto dof = rates.extent(0);

	// Both reactants are in but not the product
	if (isInSub[_reactants[0]] and isInSub[_reactants[1]]) {
		// Code not setup to deal with this
	}
	// Both reactants are out but product is in
	else if (not isInSub[_reactants[0]] and not isInSub[_reactants[1]]) {
		// Compute the flux for the 0th order moments
		double f = this->_coefs(0, 0, 0, 0) * cR1 * cR2;
		for (auto i : speciesRangeNoI) {
			f += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()] * cR2;
			f += this->_coefs(0, i() + 1, 0, 0) * cR1 * cmR2[i()];
			for (auto j : speciesRangeNoI) {
				f += this->_coefs(i() + 1, j() + 1, 0, 0) * cmR1[i()] *
					cmR2[j()];
			}
		}
		f *= this->_rate(gridIndex);

		IndexType p = 0;
		for (auto prodId : _products) {
			if (prodId == invalidIndex) {
				continue;
			}

			if (isInSub[prodId]) {
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 1 + p, 0, 0)),
					f / _productVolumes[p]);
			}
			p++;
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					f = this->_coefs(0, 0, p + 2, k() + 1) * cR1 * cR2;
					for (auto i : speciesRangeNoI) {
						f += this->_coefs(i() + 1, 0, p + 2, k() + 1) *
							cmR1[i()] * cR2;
						f += this->_coefs(0, i() + 1, p + 2, k() + 1) * cR1 *
							cmR2[i()];
						for (auto j : speciesRangeNoI) {
							f +=
								this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) *
								cmR1[i()] * cmR2[j()];
						}
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 1 + p, 1 + k(), 0)),
						f / _productVolumes[p]);
				}
			}
		}
	}
	// Only the first reactant is in not the second one
	else if (isInSub[_reactants[0]]) {
		// Compute the flux for the 0th order moments
		double f = this->_coefs(0, 0, 0, 0) * cR2;
		for (auto i : speciesRangeNoI) {
			f += this->_coefs(0, i() + 1, 0, 0) * cmR2[i()];
		}
		f *= this->_rate(gridIndex);

		// First for the first reactant
		Kokkos::atomic_sub(&rates(this->_rateEntries(subId, 0, 0, 0)),
			f / _reactantVolumes[0]);
		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}
			if (isInSub[prodId])
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 1 + p, 0, 0)),
					f / _productVolumes[p]);
		}

		// 1st moment contribution
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] == invalidIndex)
				continue;
			f = this->_coefs(i() + 1, 0, 0, 0) * cR2;
			for (auto j : speciesRangeNoI) {
				f += this->_coefs(i() + 1, j() + 1, 0, 0) * cmR2[j()];
			}
			f *= this->_rate(gridIndex);

			// First for the first reactant
			Kokkos::atomic_sub(&rates(this->_rateEntries(subId, 0, 0, 1 + i())),
				f / _reactantVolumes[0]);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (isInSub[prodId])
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 1 + p, 0, 1 + i())),
						f / _productVolumes[p]);
			}
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// First for the first reactant
			if (_reactantMomentIds[0][k()] != invalidIndex) {
				f = this->_coefs(0, 0, 0, k() + 1) * cR2;
				for (auto i : speciesRangeNoI) {
					f += this->_coefs(0, i() + 1, 0, k() + 1) * cmR2[i()];
				}
				f *= this->_rate(gridIndex);
				Kokkos::atomic_sub(
					&rates(this->_rateEntries(subId, 0, 1 + k(), 0)),
					f / _reactantVolumes[0]);

				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[0][i()] == invalidIndex)
						continue;
					f = this->_coefs(i() + 1, 0, 0, k() + 1) * cR2;
					for (auto j : speciesRangeNoI) {
						f += this->_coefs(i() + 1, j() + 1, 0, k() + 1) *
							cmR2[j()];
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_sub(
						&rates(this->_rateEntries(subId, 0, 1 + k(), 1 + i())),
						f / _reactantVolumes[0]);
				}
			}

			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					f = this->_coefs(0, 0, p + 2, k() + 1) * cR2;
					for (auto i : speciesRangeNoI) {
						f += this->_coefs(0, i() + 1, p + 2, k() + 1) *
							cmR2[i()];
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 1 + p, 1 + k(), 0)),
						f / _productVolumes[p]);

					for (auto i : speciesRangeNoI) {
						if (_reactantMomentIds[0][i()] == invalidIndex)
							continue;
						f = this->_coefs(i() + 1, 0, p + 2, k() + 1) * cR2;
						for (auto j : speciesRangeNoI) {
							f +=
								this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) *
								cmR2[j()];
						}
						f *= this->_rate(gridIndex);
						Kokkos::atomic_add(&rates(this->_rateEntries(
											   subId, 1 + p, 1 + k(), 1 + i())),
							f / _productVolumes[p]);
					}
				}
			}
		}
	}
	// Last case, only the second product is in
	else {
		// Compute the flux for the 0th order moments
		double f = this->_coefs(0, 0, 0, 0) * cR1;
		for (auto i : speciesRangeNoI) {
			f += this->_coefs(i() + 1, 0, 0, 0) * cmR1[i()];
		}
		f *= this->_rate(gridIndex);

		// First for the reactant
		Kokkos::atomic_sub(&rates(this->_rateEntries(subId, 0, 0, 0)),
			f / _reactantVolumes[1]);
		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}
			if (isInSub[prodId]) {
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 1 + p, 0, 0)),
					f / _productVolumes[p]);
			}
		}

		// Compute the flux for the 0th order moments, moment contribution
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] == invalidIndex)
				continue;
			f = this->_coefs(0, i() + 1, 0, 0) * cR1;
			for (auto j : speciesRangeNoI) {
				f += this->_coefs(i() + 1, j() + 1, 0, 0) * cmR1[i()];
			}
			f *= this->_rate(gridIndex);

			// First for the reactant
			Kokkos::atomic_sub(&rates(this->_rateEntries(subId, 0, 0, 1 + i())),
				f / _reactantVolumes[1]);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (isInSub[prodId])
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 1 + p, 0, 1 + i())),
						f / _productVolumes[p]);
			}
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the second reactant
			if (_reactantMomentIds[1][k()] != invalidIndex) {
				f = this->_coefs(0, 0, 1, k() + 1) * cR1;
				for (auto i : speciesRangeNoI) {
					f += this->_coefs(i() + 1, 0, 1, k() + 1) * cmR1[i()];
				}
				f *= this->_rate(gridIndex);
				Kokkos::atomic_sub(
					&rates(this->_rateEntries(subId, 0, 1 + k(), 0)),
					f / _reactantVolumes[1]);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[1][i()] == invalidIndex)
						continue;
					f = this->_coefs(0, i() + 1, 1, k() + 1) * cR1;
					for (auto j : speciesRangeNoI) {
						f += this->_coefs(j() + 1, i() + 1, 1, k() + 1) *
							cmR1[j()];
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_sub(
						&rates(this->_rateEntries(subId, 0, 1 + k(), 1 + i())),
						f / _reactantVolumes[1]);
				}
			}

			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					f = this->_coefs(0, 0, p + 2, k() + 1) * cR1;
					for (auto i : speciesRangeNoI) {
						f += this->_coefs(i() + 1, 0, p + 2, k() + 1) *
							cmR1[i()];
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 1 + p, 1 + k(), 0)),
						f / _productVolumes[p]);

					// 1st moment contribution
					for (auto i : speciesRangeNoI) {
						if (_reactantMomentIds[1][i()] == invalidIndex)
							continue;
						f = this->_coefs(0, i() + 1, p + 2, k() + 1) * cR1;
						for (auto j : speciesRangeNoI) {
							f +=
								this->_coefs(j() + 1, i() + 1, p + 2, k() + 1) *
								cmR1[j()];
						}
						f *= this->_rate(gridIndex);
						Kokkos::atomic_add(&rates(this->_rateEntries(
											   subId, 1 + p, 1 + k(), 1 + i())),
							f / _productVolumes[p]);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::getConstantConnectivities(
	ConnectivitiesView conns, BelongingView isInSub, OwnedSubMapView backMap)
{
	// Check products
	bool productInSub = false;
	AmountType nProd = 0;
	for (auto prodId : _products) {
		if (prodId == invalidIndex) {
			continue;
		}
		nProd++;
		if (isInSub[prodId])
			productInSub = true;
	}
	// Only consider specific cases
	if (not isInSub[_reactants[0]] and not isInSub[_reactants[1]]) {
		if (nProd == 0)
			return;
		if (nProd > 0 && not productInSub)
			return;
	}
	if (isInSub[_reactants[0]] and isInSub[_reactants[1]]) {
		if (nProd == 0)
			return;
		if (nProd > 0 && productInSub)
			return;
	}

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto dof = conns.extent(0);

	// Both reactants are in but not the product
	if (isInSub[_reactants[0]] and isInSub[_reactants[1]]) {
		// Code not setup to deal with this
	}
	// Both reactants are out but product is in
	else if (not isInSub[_reactants[0]] and not isInSub[_reactants[1]]) {
		IndexType p = 0;
		for (auto prodId : _products) {
			if (prodId == invalidIndex) {
				continue;
			}

			if (isInSub[prodId])
				conns(backMap(prodId), dof) = true;
			p++;
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					conns(backMap(_productMomentIds[p][k()]), dof) = true;
				}
			}
		}
	}
	// Only the first reactant is in not the second one
	else if (isInSub[_reactants[0]]) {
		// First for the first reactant
		conns(backMap(_reactants[0]), backMap(_reactants[0])) = true;
		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}
			if (isInSub[prodId])
				conns(backMap(prodId), backMap(_reactants[0])) = true;
		}

		// 1st moment contribution
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] == invalidIndex)
				continue;

			// First for the first reactant
			conns(backMap(_reactants[0]), backMap(_reactantMomentIds[0][i()])) =
				true;
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (isInSub[prodId])
					conns(backMap(prodId),
						backMap(_reactantMomentIds[0][i()])) = true;
			}
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// First for the first reactant
			if (_reactantMomentIds[0][k()] != invalidIndex) {
				conns(backMap(_reactantMomentIds[0][k()]),
					backMap(_reactants[0])) = true;

				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[0][i()] == invalidIndex)
						continue;
					conns(backMap(_reactantMomentIds[0][k()]),
						backMap(_reactantMomentIds[0][i()])) = true;
				}
			}

			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					conns(backMap(_productMomentIds[p][k()]),
						backMap(_reactants[0])) = true;

					for (auto i : speciesRangeNoI) {
						if (_reactantMomentIds[0][i()] == invalidIndex)
							continue;
						conns(backMap(_productMomentIds[p][k()]),
							backMap(_reactantMomentIds[0][i()])) = true;
					}
				}
			}
		}
	}
	// Last case, only the second product is in
	else {
		// First for the reactant
		conns(backMap(_reactants[1]), backMap(_reactants[1])) = true;
		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}
			if (isInSub[prodId])
				conns(backMap(prodId), backMap(_reactants[1])) = true;
		}

		// Compute the flux for the 0th order moments, moment contribution
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] == invalidIndex)
				continue;

			// First for the reactant
			conns(backMap(_reactants[1]), backMap(_reactantMomentIds[1][i()])) =
				true;
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (isInSub[prodId])
					conns(backMap(prodId),
						backMap(_reactantMomentIds[1][i()])) = true;
			}
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the second reactant
			if (_reactantMomentIds[1][k()] != invalidIndex) {
				conns(backMap(_reactantMomentIds[1][k()]),
					backMap(_reactants[1])) = true;

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[1][i()] == invalidIndex)
						continue;
					conns(backMap(_reactantMomentIds[1][k()]),
						backMap(_reactantMomentIds[1][i()])) = true;
				}
			}

			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					conns(backMap(_productMomentIds[p][k()]),
						backMap(_reactants[1])) = true;

					// 1st moment contribution
					for (auto i : speciesRangeNoI) {
						if (_reactantMomentIds[1][i()] == invalidIndex)
							continue;
						conns(backMap(_productMomentIds[p][k()]),
							backMap(_reactantMomentIds[1][i()])) = true;
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
ProductionReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// Check if our cluster is on the left side of this reaction
	if (clusterId == _reactants[0]) {
		return this->_rate(gridIndex) * concentrations[_reactants[1]] *
			this->_coefs(0, 0, 0, 0);
	}
	if (clusterId == _reactants[1]) {
		return this->_rate(gridIndex) * concentrations[_reactants[0]] *
			this->_coefs(0, 0, 0, 0);
	}

	// This cluster is not part of the reaction
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::mapJacobianEntries(
	Connectivity connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	_connEntries[0][0][0][0] = connectivity(_reactants[0], _reactants[0]);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[0][i()] != invalidIndex) {
			_connEntries[0][0][0][1 + i()] =
				connectivity(_reactants[0], _reactantMomentIds[0][i()]);
			_connEntries[0][1 + i()][0][0] =
				connectivity(_reactantMomentIds[0][i()], _reactants[0]);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[0][j()] != invalidIndex) {
					_connEntries[0][1 + i()][0][1 + j()] = connectivity(
						_reactantMomentIds[0][i()], _reactantMomentIds[0][j()]);
				}
			}
		}
	}

	_connEntries[1][0][0][0] = connectivity(_reactants[1], _reactants[0]);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[0][i()] != invalidIndex) {
			_connEntries[1][0][0][1 + i()] =
				connectivity(_reactants[1], _reactantMomentIds[0][i()]);
		}
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			_connEntries[1][1 + i()][0][0] =
				connectivity(_reactantMomentIds[1][i()], _reactants[0]);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[0][j()] != invalidIndex) {
					_connEntries[1][1 + i()][0][1 + j()] = connectivity(
						_reactantMomentIds[1][i()], _reactantMomentIds[0][j()]);
				}
			}
		}
	}

	_connEntries[0][0][1][0] = connectivity(_reactants[0], _reactants[1]);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			_connEntries[0][0][1][1 + i()] =
				connectivity(_reactants[0], _reactantMomentIds[1][i()]);
		}
		if (_reactantMomentIds[0][i()] != invalidIndex) {
			_connEntries[0][1 + i()][1][0] =
				connectivity(_reactantMomentIds[0][i()], _reactants[1]);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[1][j()] != invalidIndex) {
					_connEntries[0][1 + i()][1][1 + j()] = connectivity(
						_reactantMomentIds[0][i()], _reactantMomentIds[1][j()]);
				}
			}
		}
	}

	_connEntries[1][0][1][0] = connectivity(_reactants[1], _reactants[1]);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[1][i()] != invalidIndex) {
			_connEntries[1][0][1][1 + i()] =
				connectivity(_reactants[1], _reactantMomentIds[1][i()]);
			_connEntries[1][1 + i()][1][0] =
				connectivity(_reactantMomentIds[1][i()], _reactants[1]);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[1][j()] != invalidIndex) {
					_connEntries[1][1 + i()][1][1 + j()] = connectivity(
						_reactantMomentIds[1][i()], _reactantMomentIds[1][j()]);
				}
			}
		}
	}

	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}

		_connEntries[2 + p][0][0][0] = connectivity(prodId, _reactants[0]);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] != invalidIndex) {
				_connEntries[2 + p][0][0][1 + i()] =
					connectivity(prodId, _reactantMomentIds[0][i()]);
			}
			if (_productMomentIds[p][i()] != invalidIndex) {
				_connEntries[2 + p][1 + i()][0][0] =
					connectivity(_productMomentIds[p][i()], _reactants[0]);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[0][j()] != invalidIndex) {
						_connEntries[2 + p][1 + i()][0][1 + j()] =
							connectivity(_productMomentIds[p][i()],
								_reactantMomentIds[0][j()]);
					}
				}
			}
		}

		_connEntries[2 + p][0][1][0] = connectivity(prodId, _reactants[1]);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] != invalidIndex) {
				_connEntries[2 + p][0][1][1 + i()] =
					connectivity(prodId, _reactantMomentIds[1][i()]);
			}
			if (_productMomentIds[p][i()] != invalidIndex) {
				_connEntries[2 + p][1 + i()][1][0] =
					connectivity(_productMomentIds[p][i()], _reactants[1]);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[1][j()] != invalidIndex) {
						_connEntries[2 + p][1 + i()][1][1 + j()] =
							connectivity(_productMomentIds[p][i()],
								_reactantMomentIds[1][j()]);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::mapRateEntries(
	ConnectivitiesPairView connectivityRow,
	ConnectivitiesPairView connectivityEntries, BelongingView isInSub,
	OwnedSubMapView backMap, IndexType subId)
{
	// Check products
	bool productInSub = false;
	AmountType nProd = 0;
	for (auto prodId : _products) {
		if (prodId == invalidIndex) {
			continue;
		}
		nProd++;
		if (isInSub[prodId])
			productInSub = true;
	}
	// Only consider specific cases
	if (not isInSub[_reactants[0]] and not isInSub[_reactants[1]]) {
		if (nProd == 0)
			return;
		if (nProd > 0 && not productInSub)
			return;
	}
	if (isInSub[_reactants[0]] and isInSub[_reactants[1]]) {
		if (nProd == 0)
			return;
		if (nProd > 0 && productInSub)
			return;
	}

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto dof = connectivityRow.extent(0) - 1;

	// Both reactants are in but not the product
	if (isInSub[_reactants[0]] and isInSub[_reactants[1]]) {
		// Code not setup to deal with this
	}
	// Both reactants are out but product is in
	else if (not isInSub[_reactants[0]] and not isInSub[_reactants[1]]) {
		IndexType p = 0;
		for (auto prodId : _products) {
			if (prodId == invalidIndex) {
				continue;
			}

			if (isInSub[prodId]) {
				this->_rateEntries(subId, 1 + p, 0, 0) = this->getPosition(
					backMap(prodId), dof, connectivityRow, connectivityEntries);
			}
			p++;
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					this->_rateEntries(subId, 1 + p, 1 + k(), 0) =
						this->getPosition(backMap(_productMomentIds[p][k()]),
							dof, connectivityRow, connectivityEntries);
				}
			}
		}
	}
	// Only the first reactant is in not the second one
	else if (isInSub[_reactants[0]]) {
		// First for the first reactant
		this->_rateEntries(subId, 0, 0, 0) =
			this->getPosition(backMap(_reactants[0]), backMap(_reactants[0]),
				connectivityRow, connectivityEntries);
		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}
			if (isInSub[prodId])
				this->_rateEntries(subId, 1 + p, 0, 0) =
					this->getPosition(backMap(prodId), backMap(_reactants[0]),
						connectivityRow, connectivityEntries);
		}

		// 1st moment contribution
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[0][i()] == invalidIndex)
				continue;
			// First for the first reactant
			this->_rateEntries(subId, 0, 0, 1 + i()) = this->getPosition(
				backMap(_reactants[0]), backMap(_reactantMomentIds[0][i()]),
				connectivityRow, connectivityEntries);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (isInSub[prodId])
					this->_rateEntries(subId, 1 + p, 0, 1 + i()) =
						this->getPosition(backMap(prodId),
							backMap(_reactantMomentIds[0][i()]),
							connectivityRow, connectivityEntries);
			}
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the first reactant
			if (_reactantMomentIds[0][k()] != invalidIndex) {
				this->_rateEntries(subId, 0, 1 + k(), 0) = this->getPosition(
					backMap(_reactantMomentIds[0][k()]), backMap(_reactants[0]),
					connectivityRow, connectivityEntries);
				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[0][i()] == invalidIndex)
						continue;
					this->_rateEntries(subId, 0, 1 + k(), 1 + i()) =
						this->getPosition(backMap(_reactantMomentIds[0][k()]),
							backMap(_reactantMomentIds[0][i()]),
							connectivityRow, connectivityEntries);
				}
			}

			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					this->_rateEntries(subId, 1 + p, 1 + k(), 0) =
						this->getPosition(backMap(_productMomentIds[p][k()]),
							backMap(_reactants[0]), connectivityRow,
							connectivityEntries);
					for (auto i : speciesRangeNoI) {
						if (_reactantMomentIds[0][i()] == invalidIndex)
							continue;
						this->_rateEntries(subId, 1 + p, 1 + k(), 1 + i()) =
							this->getPosition(
								backMap(_productMomentIds[p][k()]),
								backMap(_reactantMomentIds[0][i()]),
								connectivityRow, connectivityEntries);
					}
				}
			}
		}
	}
	// Last case, only the second product is in
	else {
		// First for the reactant
		this->_rateEntries(subId, 0, 0, 0) =
			this->getPosition(backMap(_reactants[1]), backMap(_reactants[1]),
				connectivityRow, connectivityEntries);
		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}
			if (isInSub[prodId]) {
				this->_rateEntries(subId, 1 + p, 0, 0) =
					this->getPosition(backMap(prodId), backMap(_reactants[1]),
						connectivityRow, connectivityEntries);
			}
		}

		// Compute the flux for the 0th order moments, moment contribution
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[1][i()] == invalidIndex)
				continue;
			// First for the reactant
			this->_rateEntries(subId, 0, 0, 1 + i()) = this->getPosition(
				backMap(_reactants[1]), backMap(_reactantMomentIds[1][i()]),
				connectivityRow, connectivityEntries);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (isInSub[prodId])
					this->_rateEntries(subId, 1 + p, 0, 1 + i()) =
						this->getPosition(backMap(prodId),
							backMap(_reactantMomentIds[1][i()]),
							connectivityRow, connectivityEntries);
			}
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the second reactant
			if (_reactantMomentIds[1][k()] != invalidIndex) {
				this->_rateEntries(subId, 0, 1 + k(), 0) = this->getPosition(
					backMap(_reactantMomentIds[1][k()]), backMap(_reactants[1]),
					connectivityRow, connectivityEntries);
				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[1][i()] == invalidIndex)
						continue;
					this->_rateEntries(subId, 0, 1 + k(), 1 + i()) =
						this->getPosition(backMap(_reactantMomentIds[1][k()]),
							backMap(_reactantMomentIds[1][i()]),
							connectivityRow, connectivityEntries);
				}
			}

			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				if (not isInSub[prodId])
					continue;

				if (_productMomentIds[p][k()] != invalidIndex) {
					this->_rateEntries(subId, 1 + p, 1 + k(), 0) =
						this->getPosition(backMap(_productMomentIds[p][k()]),
							backMap(_reactants[1]), connectivityRow,
							connectivityEntries);
					// 1st moment contribution
					for (auto i : speciesRangeNoI) {
						if (_reactantMomentIds[1][i()] == invalidIndex)
							continue;
						this->_rateEntries(subId, 1 + p, 1 + k(), 1 + i()) =
							this->getPosition(
								backMap(_productMomentIds[p][k()]),
								backMap(_reactantMomentIds[1][i()]),
								connectivityRow, connectivityEntries);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION
DissociationReaction<TNetwork, TDerived>::DissociationReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1,
	IndexType cluster2) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_products({cluster1, cluster2})
{
	this->copyMomentIds(_reactant, _reactantMomentIds);
	for (auto i : {0, 1}) {
		this->copyMomentIds(_products[i], _productMomentIds[i]);
	}

	// static
	const auto dummyRegion = Region(Composition{});

	auto clReg = this->_clusterData->getCluster(_reactant).getRegion();
	auto prod1Reg = this->_clusterData->getCluster(_products[0]).getRegion();
	auto prod2Reg = this->_clusterData->getCluster(_products[1]).getRegion();

	_reactantVolume = clReg.volume();
	_productVolumes = {prod1Reg.volume(), prod2Reg.volume()};

	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_FUNCTION
DissociationReaction<TNetwork, TDerived>::DissociationReaction(
	ReactionDataRef reactionData, const ClusterData& clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	DissociationReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeCoefficients()
{
	// static
	const auto dummyRegion = Region(Composition{});

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto clReg = this->_clusterData->getCluster(_reactant).getRegion();
	auto prod1Reg = this->_clusterData->getCluster(_products[0]).getRegion();
	auto prod2Reg = this->_clusterData->getCluster(_products[1]).getRegion();
	const auto& clDisp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			clReg);
	const auto& prod1Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			prod1Reg);
	const auto& prod2Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			prod2Reg);
	auto cl2Reg = dummyRegion;

	// Initialize the reflected regions
	auto rRegions = detail::updateReflectedRegionsForCoefs<nMomentIds>(
		prod1Reg, prod2Reg, clReg, cl2Reg);
	auto clRR = rRegions[2];
	auto cl2RR = rRegions[3];
	auto pr1RR = rRegions[0];
	auto pr2RR = rRegions[1];

	auto nOverlap = this->computeOverlap(pr1RR, pr2RR, clRR, cl2RR);

	// The first coefficient is simply the overlap because it is the sum over 1
	this->_coefs(0, 0, 0, 0) = nOverlap;
	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		// First order sum
		this->_coefs(i() + 1, 0, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), clRR, cl2RR, pr2RR, pr1RR);
	}

	// First moments
	for (auto k : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[k()];
		// Reactant
		this->_coefs(0, 0, 0, k() + 1) =
			this->_coefs(k() + 1, 0, 0, 0) / clDisp[k()];

		// First product
		this->_coefs(0, 0, 1, k() + 1) = factor *
			detail::computeFirstOrderSum(k(), pr1RR, pr2RR, clRR, cl2RR) /
			prod1Disp[k()];

		// Second product
		this->_coefs(0, 0, 2, k() + 1) = factor *
			detail::computeFirstOrderSum(k(), pr2RR, pr1RR, clRR, cl2RR) /
			prod2Disp[k()];
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
						i(), clRR, cl2RR, pr2RR, pr1RR) /
					clDisp[k()];
				this->_coefs(i() + 1, 0, 1, k() + 1) = factor *
					detail::computeSecondOrderOffsetSum(
						i(), clRR, cl2RR, pr1RR, pr2RR) /
					prod1Disp[k()];
				this->_coefs(i() + 1, 0, 2, k() + 1) = factor *
					detail::computeSecondOrderOffsetSum(
						i(), clRR, cl2RR, pr2RR, pr1RR) /
					prod2Disp[k()];
			}
			else {
				this->_coefs(i() + 1, 0, 0, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(k() + 1, 0, 0, 0) / (nOverlap * clDisp[k()]);
				this->_coefs(i() + 1, 0, 1, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(0, 0, 1, k() + 1) / nOverlap;
				this->_coefs(i() + 1, 0, 2, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(0, 0, 2, k() + 1) / nOverlap;
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
DissociationReaction<TNetwork, TDerived>::computeRate(
	IndexType gridIndex, double time)
{
	double omega = this->_clusterData->atomicVolume();
	double T = this->_clusterData->temperature(gridIndex);
	constexpr double k_B = ::xolotl::core::kBoltzmann;

	double kPlus = this->asDerived()->getRateForProduction(gridIndex);
	double E_b = this->asDerived()->computeBindingEnergy(time);

	double kMinus = (1.0 / omega) * kPlus * std::exp(-E_b / (k_B * T));

	return kMinus;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

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
	// Each product connects with  the reactant
	// Product 1 with reactant
	this->addConnectivity(_products[0], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_products[0], _reactantMomentIds[i()], connectivity);
		}
		if (_productMomentIds[0][i()] != invalidIndex) {
			this->addConnectivity(
				_productMomentIds[0][i()], _reactant, connectivity);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					this->addConnectivity(_productMomentIds[0][i()],
						_reactantMomentIds[j()], connectivity);
				}
			}
		}
	}
	// Product 2 with reactant
	this->addConnectivity(_products[1], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_products[1], _reactantMomentIds[i()], connectivity);
		}
		if (_productMomentIds[1][i()] != invalidIndex) {
			this->addConnectivity(
				_productMomentIds[1][i()], _reactant, connectivity);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					this->addConnectivity(_productMomentIds[1][i()],
						_reactantMomentIds[j()], connectivity);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData->getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	auto prod1 = this->_clusterData->getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	auto prod2 = this->_clusterData->getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();

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
	// Each product connects with  the reactant
	// Product 1 with reactant
	if (_products[0] == _reactant)
		this->addConnectivity(_products[0], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[0][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (_productMomentIds[0][i()] == _reactantMomentIds[j()])
					this->addConnectivity(_productMomentIds[0][i()],
						_reactantMomentIds[j()], connectivity);
			}
		}
	}
	// Product 2 with reactant
	if (_products[1] == _reactant)
		this->addConnectivity(_products[1], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[1][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (_productMomentIds[1][i()] == _reactantMomentIds[j()])
					this->addConnectivity(_productMomentIds[1][i()],
						_reactantMomentIds[j()], connectivity);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeFlux(
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
	Kokkos::atomic_add(&fluxes[_products[0]], f / _productVolumes[0]);
	Kokkos::atomic_add(&fluxes[_products[1]], f / _productVolumes[1]);

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

		// Now the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 1, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 1, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(
				&fluxes[_productMomentIds[0][k()]], f / _productVolumes[0]);
		}

		// Finally the second product
		if (_productMomentIds[1][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 2, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 2, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(
				&fluxes[_productMomentIds[1][k()]], f / _productVolumes[1]);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computePartialDerivatives(
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
	// For the first product
	df = this->_rate(gridIndex) / _productVolumes[0];
	Kokkos::atomic_add(
		&values(_connEntries[1][0][0][0]), df * this->_coefs(0, 0, 0, 0));

	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_add(&values(_connEntries[1][0][0][1 + i()]),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the second product
	df = this->_rate(gridIndex) / _productVolumes[1];
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
		// For the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / _productVolumes[0];
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
			df = this->_rate(gridIndex) / _productVolumes[1];
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

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
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
	df = this->_rate(gridIndex) / _productVolumes[0];
	if (_products[0] == _reactant)
		Kokkos::atomic_add(
			&values(_connEntries[1][0][0][0]), df * this->_coefs(0, 0, 0, 0));

	// For the second product
	df = this->_rate(gridIndex) / _productVolumes[1];
	if (_products[1] == _reactant)
		Kokkos::atomic_add(
			&values(_connEntries[2][0][0][0]), df * this->_coefs(0, 0, 0, 0));

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
		// For the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / _productVolumes[0];
			for (auto i : speciesRangeNoI) {
				if (_productMomentIds[0][k()] == _reactantMomentIds[i()])
					Kokkos::atomic_add(
						&values(_connEntries[1][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 1, k() + 1));
			}
		}
		// For the second product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / _productVolumes[1];
			for (auto i : speciesRangeNoI) {
				if (_productMomentIds[1][k()] == _reactantMomentIds[i()])
					Kokkos::atomic_add(
						&values(_connEntries[2][1 + k()][0][1 + i()]),
						df * this->_coefs(i() + 1, 0, 2, k() + 1));
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeConstantRates(
	ConcentrationsView concentrations, RatesView rates, BelongingView isInSub,
	IndexType subId, IndexType gridIndex)
{
	// Only consider cases specific cases
	if (not isInSub[_reactant] and not isInSub[_products[0]] and
		not isInSub[_products[1]])
		return;
	if (isInSub[_reactant] and isInSub[_products[0]] and isInSub[_products[1]])
		return;

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

	// Compute the terms for the 0th order moments
	// First case where the reactant is in
	if (isInSub[_reactant]) {
		// Compute the flux for the 0th order moments
		double f = this->_coefs(0, 0, 0, 0) * this->_rate(gridIndex);
		Kokkos::atomic_sub(
			&rates(this->_rateEntries(subId, 0, 0, 0)), f / _reactantVolume);
		if (isInSub[_products[0]])
			Kokkos::atomic_add(&rates(this->_rateEntries(subId, 1, 0, 0)),
				f / _productVolumes[0]);
		if (isInSub[_products[1]])
			Kokkos::atomic_add(&rates(this->_rateEntries(subId, 2, 0, 0)),
				f / _productVolumes[1]);

		// Now the moment contribtions
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] == invalidIndex)
				continue;

			f = this->_coefs(i() + 1, 0, 0, 0) * this->_rate(gridIndex);
			Kokkos::atomic_sub(&rates(this->_rateEntries(subId, 0, 0, 1 + i())),
				f / _reactantVolume);
			if (isInSub[_products[0]])
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 1, 0, 1 + i())),
					f / _productVolumes[0]);
			if (isInSub[_products[1]])
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 2, 0, 1 + i())),
					f / _productVolumes[1]);
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// First for the reactant
			if (_reactantMomentIds[k()] != invalidIndex) {
				f = this->_coefs(0, 0, 0, k() + 1) * this->_rate(gridIndex);
				Kokkos::atomic_sub(
					&rates(this->_rateEntries(subId, 0, 1 + k(), 0)),
					f / _reactantVolume);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					f = this->_coefs(i() + 1, 0, 0, k() + 1) *
						this->_rate(gridIndex);
					Kokkos::atomic_sub(
						&rates(this->_rateEntries(subId, 0, 1 + k(), 1 + i())),
						f / _reactantVolume);
				}
			}

			// Now the first product
			if (isInSub[_products[0]] and
				_productMomentIds[0][k()] != invalidIndex) {
				f = this->_coefs(0, 0, 1, k() + 1) * this->_rate(gridIndex);
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 1, 1 + k(), 0)),
					f / _productVolumes[0]);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					f = this->_coefs(i() + 1, 0, 1, k() + 1) *
						this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 0, 1 + k(), 1 + i())),
						f / _productVolumes[0]);
				}
			}

			// Finally the second product
			if (isInSub[_products[1]] and
				_productMomentIds[1][k()] != invalidIndex) {
				f = this->_coefs(0, 0, 2, k() + 1) * this->_rate(gridIndex);
				Kokkos::atomic_add(
					&rates(this->_rateEntries(subId, 2, 1 + k(), 0)),
					f / _productVolumes[1]);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					f = this->_coefs(i() + 1, 0, 2, k() + 1) *
						this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 2, 1 + k(), 1 + i())),
						f / _productVolumes[1]);
				}
			}
		}
	}
	// Now the reactant is not in
	else {
		auto dof = rates.extent(0);
		double f = this->_coefs(0, 0, 0, 0) * cR;
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] != invalidIndex)
				f += this->_coefs(i() + 1, 0, 0, 0) * cmR[i()];
		}
		f *= this->_rate(gridIndex);

		// For the first product
		if (isInSub[_products[0]])
			Kokkos::atomic_add(&rates(this->_rateEntries(subId, 1, 0, 0)),
				f / (double)_productVolumes[0]);
		// For the second product
		if (isInSub[_products[1]])
			Kokkos::atomic_add(&rates(this->_rateEntries(subId, 2, 0, 0)),
				f / (double)_productVolumes[1]);

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the first product
			if (isInSub[_products[0]]) {
				if (_productMomentIds[0][k()] != invalidIndex) {
					f = this->_coefs(0, 0, 1, k() + 1) * cR;
					for (auto i : speciesRangeNoI) {
						f += this->_coefs(i() + 1, 0, 1, k() + 1) * cmR[i()];
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 1, 1 + k(), 0)),
						f / _productVolumes[0]);
				}
			}

			// For the second product
			if (isInSub[_products[1]]) {
				if (_productMomentIds[1][k()] != invalidIndex) {
					f = this->_coefs(0, 0, 2, k() + 1) * cR;
					for (auto i : speciesRangeNoI) {
						f += this->_coefs(i() + 1, 0, 2, k() + 1) * cmR[i()];
					}
					f *= this->_rate(gridIndex);
					Kokkos::atomic_add(
						&rates(this->_rateEntries(subId, 2, 1 + k(), 0)),
						f / _productVolumes[1]);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::getConstantConnectivities(
	ConnectivitiesView conns, BelongingView isInSub, OwnedSubMapView backMap)
{
	// Only consider cases specific cases
	if (not isInSub[_reactant] and not isInSub[_products[0]] and
		not isInSub[_products[1]])
		return;
	if (isInSub[_reactant] and isInSub[_products[0]] and isInSub[_products[1]])
		return;

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Terms for the 0th order moments
	// First case where the reactant is in
	if (isInSub[_reactant]) {
		conns(backMap(_reactant), backMap(_reactant)) = true;
		if (isInSub[_products[0]])
			conns(backMap(_products[0]), backMap(_reactant)) = true;
		if (isInSub[_products[1]])
			conns(backMap(_products[1]), backMap(_reactant)) = true;

		// Now the moment contribtions
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] == invalidIndex)
				continue;

			conns(backMap(_reactant), backMap(_reactantMomentIds[i()])) = true;
			if (isInSub[_products[0]])
				conns(backMap(_products[0]), backMap(_reactantMomentIds[i()])) =
					true;
			if (isInSub[_products[1]])
				conns(backMap(_products[1]), backMap(_reactantMomentIds[i()])) =
					true;
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// First for the reactant
			if (_reactantMomentIds[k()] != invalidIndex) {
				conns(backMap(_reactantMomentIds[k()]), backMap(_reactant)) =
					true;

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					conns(backMap(_reactantMomentIds[k()]),
						backMap(_reactantMomentIds[i()])) = true;
				}
			}

			// Now the first product
			if (isInSub[_products[0]] and
				_productMomentIds[0][k()] != invalidIndex) {
				conns(backMap(_productMomentIds[0][k()]), backMap(_reactant)) =
					true;

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					conns(backMap(_productMomentIds[0][k()]),
						backMap(_reactantMomentIds[i()])) = true;
				}
			}

			// Finally the second product
			if (isInSub[_products[1]] and
				_productMomentIds[1][k()] != invalidIndex) {
				conns(backMap(_productMomentIds[1][k()]), backMap(_reactant)) =
					true;

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					conns(backMap(_productMomentIds[1][k()]),
						backMap(_reactantMomentIds[i()])) = true;
				}
			}
		}
	}
	// Now the reactant is not in
	else {
		auto dof = conns.extent(0);

		// For the first product
		if (isInSub[_products[0]])
			conns(backMap(_products[0]), dof) = true;
		// For the second product
		if (isInSub[_products[1]])
			conns(backMap(_products[1]), dof) = true;

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the first product
			if (isInSub[_products[0]]) {
				if (_productMomentIds[0][k()] != invalidIndex) {
					conns(backMap(_productMomentIds[0][k()]), dof) = true;
				}
			}

			// For the second product
			if (isInSub[_products[1]]) {
				if (_productMomentIds[1][k()] != invalidIndex) {
					conns(backMap(_productMomentIds[1][k()]), dof) = true;
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
DissociationReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// Check if our cluster is on the left side of this reaction
	if (clusterId == _reactant) {
		return this->_rate(gridIndex) * this->_coefs(0, 0, 0, 0);
	}

	// This cluster is not part of the reaction
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::mapJacobianEntries(
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

	for (auto p : {0, 1}) {
		_connEntries[1 + p][0][0][0] = connectivity(_products[p], _reactant);
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] != invalidIndex) {
				_connEntries[1 + p][0][0][1 + i()] =
					connectivity(_products[p], _reactantMomentIds[i()]);
			}
			if (_productMomentIds[p][i()] != invalidIndex) {
				_connEntries[1 + p][1 + i()][0][0] =
					connectivity(_productMomentIds[p][i()], _reactant);
				for (auto j : speciesRangeNoI) {
					if (_reactantMomentIds[j()] != invalidIndex) {
						_connEntries[1 + p][1 + i()][0][1 + j()] = connectivity(
							_productMomentIds[p][i()], _reactantMomentIds[j()]);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::mapRateEntries(
	ConnectivitiesPairView connectivityRow,
	ConnectivitiesPairView connectivityEntries, BelongingView isInSub,
	OwnedSubMapView backMap, IndexType subId)
{
	// Only consider cases specific cases
	if (not isInSub[_reactant] and not isInSub[_products[0]] and
		not isInSub[_products[1]])
		return;
	if (isInSub[_reactant] and isInSub[_products[0]] and isInSub[_products[1]])
		return;

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the terms for the 0th order moments
	// First case where the reactant is in
	if (isInSub[_reactant]) {
		// Compute the flux for the 0th order moments
		this->_rateEntries(subId, 0, 0, 0) =
			this->getPosition(backMap(_reactant), backMap(_reactant),
				connectivityRow, connectivityEntries);
		if (isInSub[_products[0]])
			this->_rateEntries(subId, 1, 0, 0) =
				this->getPosition(backMap(_products[0]), backMap(_reactant),
					connectivityRow, connectivityEntries);
		if (isInSub[_products[1]])
			this->_rateEntries(subId, 2, 0, 0) =
				this->getPosition(backMap(_products[1]), backMap(_reactant),
					connectivityRow, connectivityEntries);

		// Now the moment contributions
		for (auto i : speciesRangeNoI) {
			if (_reactantMomentIds[i()] == invalidIndex)
				continue;

			this->_rateEntries(subId, 0, 0, 1 + i()) = this->getPosition(
				backMap(_reactant), backMap(_reactantMomentIds[i()]),
				connectivityRow, connectivityEntries);
			if (isInSub[_products[0]])
				this->_rateEntries(subId, 1, 0, 1 + i()) = this->getPosition(
					backMap(_products[0]), backMap(_reactantMomentIds[i()]),
					connectivityRow, connectivityEntries);
			if (isInSub[_products[1]])
				this->_rateEntries(subId, 2, 0, 1 + i()) = this->getPosition(
					backMap(_products[1]), backMap(_reactantMomentIds[i()]),
					connectivityRow, connectivityEntries);
		}

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// First for the reactant
			if (_reactantMomentIds[k()] != invalidIndex) {
				this->_rateEntries(subId, 0, 1 + k(), 0) = this->getPosition(
					backMap(_reactantMomentIds[k()]), backMap(_reactant),
					connectivityRow, connectivityEntries);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					this->_rateEntries(subId, 0, 1 + k(), 1 + i()) =
						this->getPosition(backMap(_reactantMomentIds[k()]),
							backMap(_reactantMomentIds[i()]), connectivityRow,
							connectivityEntries);
				}
			}

			// Now the first product
			if (isInSub[_products[0]] and
				_productMomentIds[0][k()] != invalidIndex) {
				this->_rateEntries(subId, 1, 1 + k(), 0) = this->getPosition(
					backMap(_productMomentIds[0][k()]), backMap(_reactant),
					connectivityRow, connectivityEntries);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					this->_rateEntries(subId, 1, 1 + k(), 1 + i()) =
						this->getPosition(backMap(_productMomentIds[0][k()]),
							backMap(_reactantMomentIds[i()]), connectivityRow,
							connectivityEntries);
				}
			}

			// Finally the second product
			if (isInSub[_products[1]] and
				_productMomentIds[1][k()] != invalidIndex) {
				this->_rateEntries(subId, 2, 1 + k(), 0) = this->getPosition(
					backMap(_productMomentIds[1][k()]), backMap(_reactant),
					connectivityRow, connectivityEntries);

				// 1st moment contribution
				for (auto i : speciesRangeNoI) {
					if (_reactantMomentIds[i()] == invalidIndex)
						continue;
					this->_rateEntries(subId, 2, 1 + k(), 1 + i()) =
						this->getPosition(backMap(_productMomentIds[1][k()]),
							backMap(_reactantMomentIds[i()]), connectivityRow,
							connectivityEntries);
				}
			}
		}
	}
	// Now the reactant is not in
	else {
		auto dof = connectivityRow.extent(0) - 1;

		// For the first product
		if (isInSub[_products[0]])
			this->_rateEntries(subId, 1, 0, 0) =
				this->getPosition(backMap(_products[0]), dof, connectivityRow,
					connectivityEntries);
		// For the second product
		if (isInSub[_products[1]])
			this->_rateEntries(subId, 2, 0, 0) =
				this->getPosition(backMap(_products[1]), dof, connectivityRow,
					connectivityEntries);

		// Take care of the first moments
		for (auto k : speciesRangeNoI) {
			// For the first product
			if (isInSub[_products[0]]) {
				if (_productMomentIds[0][k()] != invalidIndex) {
					this->_rateEntries(subId, 1, 1 + k(), 0) =
						this->getPosition(backMap(_productMomentIds[0][k()]),
							dof, connectivityRow, connectivityEntries);
				}
			}

			// For the second product
			if (isInSub[_products[1]]) {
				if (_productMomentIds[1][k()] != invalidIndex) {
					this->_rateEntries(subId, 2, 1 + k(), 0) =
						this->getPosition(backMap(_productMomentIds[1][k()]),
							dof, connectivityRow, connectivityEntries);
				}
			}
		}
	}
}
} // namespace network
} // namespace core
} // namespace xolotl
